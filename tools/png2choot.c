/*
 * png2choot -- PNG to CHOOT v0 encoder  (high-quality, multi-stage)
 *
 * Algorithm overview
 * ------------------
 *   Stage 1 - Initial clustering
 *     1. Load PNG via WIC -> linear RGB.
 *     2. Build weighted 5-D features: (w*nx, w*ny, R, G, B).
 *     3. k-means++ seed -> Lloyd refine -> fit Gaussian atoms.
 *
 *   Stage 2 - Iterative refinement (analysis-by-synthesis)
 *     For each refinement pass:
 *       a. Render current atoms -> float image.
 *       b. Compute per-pixel squared error vs. ground truth.
 *       c. Accumulate error per atom (each pixel blames its
 *          strongest-contributing atom).
 *       d. Split the worst N atoms: replace each with two children
 *          offset along the major eigenvector with halved variance.
 *       e. Re-run Lloyd on the enlarged atom set to settle.
 *       f. Re-fit all atoms from pixel assignments.
 *
 *   Stage 3 - Final polish
 *     After all splits, run soft colour re-fitting:
 *     for each atom, compute the mean colour of all nearby pixels
 *     weighted by the Gaussian itself (soft assignment), so the
 *     atom colour matches what the *renderer* would produce.
 *
 * This produces dramatically better quality than a single k-means pass.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "choot.h"
#include "choot_render.h"
#include "wic_png.h"

/* ------------------------------------------------------------------ */
/* helpers                                                             */
/* ------------------------------------------------------------------ */

static float clamp01(float v) {
    if (v < 0.0f) return 0.0f;
    if (v > 1.0f) return 1.0f;
    return v;
}

static float srgb_to_linear(float c) {
    if (c <= 0.04045f) return c / 12.92f;
    return powf((c + 0.055f) / 1.055f, 2.4f);
}

/* ------------------------------------------------------------------ */
/* 5-D feature point                                                   */
/* ------------------------------------------------------------------ */

typedef struct Pixel5D {
    float f[5];   /* w*nx, w*ny, r, g, b */
    float a;
} Pixel5D;

static float dist2_5d(const float a[5], const float b[5]) {
    float s = 0.0f;
    for (int i = 0; i < 5; ++i) {
        float d = a[i] - b[i];
        s += d * d;
    }
    return s;
}

/* deterministic xorshift32 */
static uint32_t xorshift32(uint32_t* state) {
    uint32_t x = *state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *state = x;
    return x;
}

/* ------------------------------------------------------------------ */
/* k-means++ seeding                                                   */
/* ------------------------------------------------------------------ */

static void kmeans_pp_seed(const Pixel5D* pixels, int n,
                           float* centers, int k, float* dist_buf) {
    uint32_t rng = (uint32_t)n ^ 0xDEADu;

    /* first centre: closest to global mean */
    {
        double mean[5] = {0};
        for (int i = 0; i < n; ++i)
            for (int d = 0; d < 5; ++d)
                mean[d] += pixels[i].f[d];
        for (int d = 0; d < 5; ++d) mean[d] /= n;

        float best = FLT_MAX;
        int best_i = 0;
        for (int i = 0; i < n; ++i) {
            float dd = 0.0f;
            for (int d = 0; d < 5; ++d) {
                float v = pixels[i].f[d] - (float)mean[d];
                dd += v * v;
            }
            if (dd < best) { best = dd; best_i = i; }
        }
        memcpy(centers, pixels[best_i].f, 5 * sizeof(float));
    }

    for (int i = 0; i < n; ++i)
        dist_buf[i] = dist2_5d(pixels[i].f, centers);

    for (int c = 1; c < k; ++c) {
        double total = 0.0;
        for (int i = 0; i < n; ++i) total += (double)dist_buf[i];
        double threshold = ((double)(xorshift32(&rng) & 0x7FFFFFFF)
                            / (double)0x7FFFFFFF) * total;
        double running = 0.0;
        int chosen = n - 1;
        for (int i = 0; i < n; ++i) {
            running += (double)dist_buf[i];
            if (running >= threshold) { chosen = i; break; }
        }
        memcpy(centers + c * 5, pixels[chosen].f, 5 * sizeof(float));
        { int i;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (i = 0; i < n; ++i) {
            float d = dist2_5d(pixels[i].f, centers + c * 5);
            if (d < dist_buf[i]) dist_buf[i] = d;
        } }
    }
}

/* ------------------------------------------------------------------ */
/* Lloyd iterations                                                    */
/* ------------------------------------------------------------------ */

static void kmeans_lloyd(const Pixel5D* pixels, int n,
                         float* centers, int k,
                         int* assignments, int max_iters) {
    double* accum  = (double*)calloc((size_t)k * 5, sizeof(double));
    int*    counts = (int*)calloc((size_t)k, sizeof(int));

    for (int iter = 0; iter < max_iters; ++iter) {
        int changed = 0;
        { int i;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(|:changed)
#endif
        for (i = 0; i < n; ++i) {
            float best = FLT_MAX;
            int best_c = 0;
            for (int c = 0; c < k; ++c) {
                float d = dist2_5d(pixels[i].f, centers + c * 5);
                if (d < best) { best = d; best_c = c; }
            }
            if (assignments[i] != best_c) {
                assignments[i] = best_c;
                changed = 1;
            }
        } }
        if (!changed && iter > 0) break;

        memset(accum,  0, (size_t)k * 5 * sizeof(double));
        memset(counts, 0, (size_t)k * sizeof(int));
        for (int i = 0; i < n; ++i) {
            int c = assignments[i];
            counts[c]++;
            for (int d = 0; d < 5; ++d)
                accum[c * 5 + d] += (double)pixels[i].f[d];
        }

        /* Reseed empty clusters: move to the farthest pixel from its
         * current centre.  This prevents dead atoms entirely. */
        {
            /* first, compute per-pixel distance to its assigned centre */
            float* pdist = (float*)malloc((size_t)n * sizeof(float));
            if (pdist) {
                { int i;
                for (i = 0; i < n; ++i) {
                    pdist[i] = dist2_5d(pixels[i].f, centers + assignments[i] * 5);
                } }

                for (int c = 0; c < k; ++c) {
                    if (counts[c] > 0) continue;
                    /* find the pixel farthest from its centre */
                    float max_d = -1.0f;
                    int max_i = 0;
                    { int i;
                    for (i = 0; i < n; ++i) {
                        if (pdist[i] > max_d) { max_d = pdist[i]; max_i = i; }
                    } }
                    /* steal this pixel as the new centre */
                    memcpy(centers + c * 5, pixels[max_i].f, 5 * sizeof(float));
                    for (int d = 0; d < 5; ++d)
                        accum[c * 5 + d] = (double)pixels[max_i].f[d];
                    counts[c] = 1;
                    /* reassign that pixel and suppress its distance */
                    assignments[max_i] = c;
                    pdist[max_i] = 0.0f;
                }
                free(pdist);
            }
        }

        for (int c = 0; c < k; ++c) {
            if (counts[c] > 0) {
                for (int d = 0; d < 5; ++d)
                    centers[c * 5 + d] = (float)(accum[c * 5 + d] / counts[c]);
            }
        }
    }

    free(accum);
    free(counts);
}

/* ------------------------------------------------------------------ */
/* Fit atoms from cluster assignments (improved)                       */
/* ------------------------------------------------------------------ */

static void fit_atoms(const Pixel5D* pixels, int n,
                      const int* assignments, int k,
                      float spatial_weight,
                      float sharpness,
                      uint32_t img_w, uint32_t img_h,
                      ChootAtom* atoms) {
    double* sx  = (double*)calloc((size_t)k, sizeof(double));
    double* sy  = (double*)calloc((size_t)k, sizeof(double));
    double* sr  = (double*)calloc((size_t)k, sizeof(double));
    double* sg  = (double*)calloc((size_t)k, sizeof(double));
    double* sb  = (double*)calloc((size_t)k, sizeof(double));
    double* sa  = (double*)calloc((size_t)k, sizeof(double));
    int*    cnt = (int*)calloc((size_t)k, sizeof(int));

    float inv_sw = 1.0f / spatial_weight;

    /* Note: scatter-add by assignment - not trivially parallelizable
     * without per-thread accumulators.  We use them here. */
#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
#else
    int nthreads = 1;
#endif
    double* t_sx  = (double*)calloc((size_t)nthreads * k, sizeof(double));
    double* t_sy  = (double*)calloc((size_t)nthreads * k, sizeof(double));
    double* t_sr  = (double*)calloc((size_t)nthreads * k, sizeof(double));
    double* t_sg  = (double*)calloc((size_t)nthreads * k, sizeof(double));
    double* t_sb  = (double*)calloc((size_t)nthreads * k, sizeof(double));
    double* t_sa  = (double*)calloc((size_t)nthreads * k, sizeof(double));
    int*    t_cnt = (int*)calloc((size_t)nthreads * k, sizeof(int));

    { int i;
#ifdef _OPENMP
#pragma omp parallel
{
    int tid = omp_get_thread_num();
#pragma omp for schedule(static)
#else
    int tid = 0;
    {
#endif
    for (i = 0; i < n; ++i) {
        int c = assignments[i];
        float nx = pixels[i].f[0] * inv_sw;
        float ny = pixels[i].f[1] * inv_sw;
        t_sx[tid * k + c] += nx;  t_sy[tid * k + c] += ny;
        t_sr[tid * k + c] += pixels[i].f[2];
        t_sg[tid * k + c] += pixels[i].f[3];
        t_sb[tid * k + c] += pixels[i].f[4];
        t_sa[tid * k + c] += pixels[i].a;
        t_cnt[tid * k + c]++;
    }
    } }
    /* reduce */
    for (int t = 0; t < nthreads; ++t) {
        for (int c = 0; c < k; ++c) {
            sx[c]  += t_sx[t * k + c];
            sy[c]  += t_sy[t * k + c];
            sr[c]  += t_sr[t * k + c];
            sg[c]  += t_sg[t * k + c];
            sb[c]  += t_sb[t * k + c];
            sa[c]  += t_sa[t * k + c];
            cnt[c] += t_cnt[t * k + c];
        }
    }
    free(t_sx); free(t_sy); free(t_sr); free(t_sg); free(t_sb); free(t_sa); free(t_cnt);

    for (int c = 0; c < k; ++c) {
        if (cnt[c] == 0) continue;
        sx[c] /= cnt[c];  sy[c] /= cnt[c];
    }

    /* covariance */
    double* cov_xx = (double*)calloc((size_t)k, sizeof(double));
    double* cov_xy = (double*)calloc((size_t)k, sizeof(double));
    double* cov_yy = (double*)calloc((size_t)k, sizeof(double));

    for (int i = 0; i < n; ++i) {
        int c = assignments[i];
        if (cnt[c] == 0) continue;
        double dx = (double)(pixels[i].f[0] * inv_sw) - sx[c];
        double dy = (double)(pixels[i].f[1] * inv_sw) - sy[c];
        cov_xx[c] += dx * dx;
        cov_xy[c] += dx * dy;
        cov_yy[c] += dy * dy;
    }

    /*
     * Variance floor: 1/4 pixel in normalised coords.
     * Prevents degenerate zero-area splats for tiny clusters.
     */
    float min_var_x = 1.0f / (float)img_w;
    float min_var_y = 1.0f / (float)img_h;
    float min_sxx = min_var_x * min_var_x;
    float min_syy = min_var_y * min_var_y;

    for (int c = 0; c < k; ++c) {
        ChootAtom* a = &atoms[c];
        memset(a, 0, sizeof(*a));

        if (cnt[c] == 0) {
            a->alpha = 0.0f;
            a->sxx = min_sxx;
            a->syy = min_syy;
            continue;
        }

        double inv_n = 1.0 / cnt[c];

        a->x = clamp01((float)sx[c]);
        a->y = clamp01((float)sy[c]);

        float mean_r = (float)(sr[c] * inv_n);
        float mean_g = (float)(sg[c] * inv_n);
        float mean_b = (float)(sb[c] * inv_n);
        float mean_a = (float)(sa[c] * inv_n);

        float luma, co, cg;
        choot_rgb_to_ycocg(mean_r, mean_g, mean_b, &luma, &co, &cg);
        a->Y     = luma;
        a->Co    = co;
        a->Cg    = cg;
        /* For mostly-opaque images, set alpha to 1.0.
         * The renderer normalises by sum(w), so alpha controls
         * relative strength, not transparency.  Using 1.0 gives
         * the strongest possible signal per atom. */
        a->alpha = (mean_a > 0.5f) ? 1.0f : clamp01(mean_a * 2.0f);

        /*
         * Scale covariance. We use the full sample covariance
         * (scale 1.0) to ensure good spatial coverage. Atoms need
         * to overlap for the weighted average to produce smooth
         * results without black gaps.
         */
        float scale = (sharpness > 0.0f) ? (1.0f / sharpness) : 1.0f;
        if (scale < 0.35f) scale = 0.35f;
        if (scale > 1.5f)  scale = 1.5f;
        float sxx = (float)(cov_xx[c] * inv_n) * scale;
        float sxy = (float)(cov_xy[c] * inv_n) * scale;
        float syy = (float)(cov_yy[c] * inv_n) * scale;

        if (sxx < min_sxx) sxx = min_sxx;
        if (syy < min_syy) syy = min_syy;

        float det = sxx * syy - sxy * sxy;
        if (det < min_sxx * min_syy * 0.01f) {
            float lim = sqrtf(sxx * syy) * 0.99f;
            sxy = (sxy > 0.0f) ? lim : -lim;
        }

        a->sxx = sxx;
        a->sxy = sxy;
        a->syy = syy;
        a->flags = 0;
    }

    free(sx); free(sy);
    free(sr); free(sg); free(sb); free(sa);
    free(cnt);
    free(cov_xx); free(cov_xy); free(cov_yy);
}

/* ------------------------------------------------------------------ */
/* Soft colour re-fitting (analysis-by-synthesis)                      */
/*                                                                     */
/* For each atom, accumulate the Gaussian-weighted mean colour of all  */
/* nearby pixels.  This makes the atom colour match what the renderer  */
/* would actually produce, instead of the hard-cluster average.        */
/* ------------------------------------------------------------------ */

static void soft_refit_colours(const Pixel5D* pixels, int n,
                               uint32_t img_w, uint32_t img_h,
                               ChootAtom* atoms, int k,
                               float spatial_weight,
                               int sample_stride) {
    (void)n;
    float inv_sw = 1.0f / spatial_weight;
    (void)inv_sw;
    if (sample_stride < 1) sample_stride = 1;
    double sample_w = (double)sample_stride * (double)sample_stride;

    { int c;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16)
#endif
    for (c = 0; c < k; ++c) {
        if (atoms[c].alpha <= 0.0f) continue;

        ChootAtom* a = &atoms[c];
        float det = a->sxx * a->syy - a->sxy * a->sxy;
        if (det <= 0.0f) continue;
        float inv_det = 1.0f / det;
        float inv_xx =  a->syy * inv_det;
        float inv_xy = -a->sxy * inv_det;
        float inv_yy =  a->sxx * inv_det;

        double wr = 0, wg = 0, wb = 0, wa_sum = 0;

        for (uint32_t py = 0; py < img_h; py += (uint32_t)sample_stride) {
            for (uint32_t px = 0; px < img_w; px += (uint32_t)sample_stride) {
                int i = (int)(py * img_w + px);
                /* pixel position in normalised coords */
                float pxn = pixels[i].f[0] / spatial_weight;
                float pyn = pixels[i].f[1] / spatial_weight;
                float dx = pxn - a->x;
                float dy = pyn - a->y;
                float quad = dx * (inv_xx * dx + inv_xy * dy)
                           + dy * (inv_xy * dx + inv_yy * dy);
                if (quad > 12.0f) continue;  /* skip negligible contributions */
                float g = expf(-0.5f * quad);
                float w = a->alpha * g;
                if (w < 1e-6f) continue;
                w *= (float)sample_w;

                wr += w * pixels[i].f[2];
                wg += w * pixels[i].f[3];
                wb += w * pixels[i].f[4];
                wa_sum += w;
            }
        }

        if (wa_sum > 1e-8) {
            float mr = (float)(wr / wa_sum);
            float mg = (float)(wg / wa_sum);
            float mb = (float)(wb / wa_sum);
            float luma, co, cg;
            choot_rgb_to_ycocg(mr, mg, mb, &luma, &co, &cg);
            a->Y  = luma;
            a->Co = co;
            a->Cg = cg;
        }
    } }
}

/* ------------------------------------------------------------------ */
/* Error-driven atom splitting                                         */
/*                                                                     */
/* 1. Render current atoms at input resolution.                        */
/* 2. Compute per-pixel MSE vs ground truth.                           */
/* 3. For each atom, accumulate error from its Voronoi cell.           */
/* 4. Split the worst atoms: offset children along major axis.         */
/* ------------------------------------------------------------------ */

typedef struct AtomError {
    int    idx;
    double error;
    int    count;
} AtomError;

static int cmp_atom_error_desc(const void* a, const void* b) {
    double ea = ((const AtomError*)a)->error;
    double eb = ((const AtomError*)b)->error;
    if (ea > eb) return -1;
    if (ea < eb) return  1;
    return 0;
}

static int split_worst_atoms(
    const Pixel5D* pixels, int n,
    uint32_t img_w, uint32_t img_h,
    ChootAtom** p_atoms, int* p_k,
    int max_split,
    int max_total_atoms,
    float spatial_weight,
    int sample_stride)
{
    (void)n;
    (void)spatial_weight;
    int k = *p_k;
    ChootAtom* atoms = *p_atoms;

    /* render current state at image resolution */
    float* rendered = (float*)malloc((size_t)img_w * img_h * 3 * sizeof(float));
    if (!rendered) return 0;

    ChootImage tmp;
    tmp.atom_count = (uint32_t)k;
    tmp.atoms = atoms;
    choot_render_linear_rgb(&tmp, (int)img_w, (int)img_h, rendered);

    /* per-atom error accumulation */
    AtomError* aerr = (AtomError*)calloc((size_t)k, sizeof(AtomError));

    /* for each pixel, find its strongest atom and accumulate error */
    if (sample_stride < 1) sample_stride = 1;
    double sample_w = (double)sample_stride * (double)sample_stride;

    for (int py = 0; py < (int)img_h; py += sample_stride) {
        for (int px = 0; px < (int)img_w; px += sample_stride) {
            int idx = py * (int)img_w + px;
            float pnx = ((float)px + 0.5f) / (float)img_w;
            float pny = ((float)py + 0.5f) / (float)img_h;

            /* find atom with strongest contribution to this pixel */
            float max_w = 0.0f;
            int max_atom = 0;
            for (int c = 0; c < k; ++c) {
                if (atoms[c].alpha <= 0.0f) continue;
                float dx = pnx - atoms[c].x;
                float dy = pny - atoms[c].y;
                float det = atoms[c].sxx * atoms[c].syy - atoms[c].sxy * atoms[c].sxy;
                if (det <= 0.0f) continue;
                float inv_det = 1.0f / det;
                float quad = dx * (atoms[c].syy * inv_det * dx + (-atoms[c].sxy * inv_det) * dy)
                           + dy * ((-atoms[c].sxy * inv_det) * dx + atoms[c].sxx * inv_det * dy);
                if (quad > 16.0f) continue;
                float g = expf(-0.5f * quad);
                float w = atoms[c].alpha * g;
                if (w > max_w) { max_w = w; max_atom = c; }
            }

            /* pixel error */
            float gt_r = pixels[idx].f[2];
            float gt_g = pixels[idx].f[3];
            float gt_b = pixels[idx].f[4];
            float re_r = rendered[idx * 3 + 0];
            float re_g = rendered[idx * 3 + 1];
            float re_b = rendered[idx * 3 + 2];
            float dr = gt_r - re_r;
            float dg = gt_g - re_g;
            float db = gt_b - re_b;
            float err = dr * dr + dg * dg + db * db;
            err *= (float)sample_w;

            aerr[max_atom].error += err;
            aerr[max_atom].count++;
        }
    }

    for (int c = 0; c < k; ++c) aerr[c].idx = c;
    qsort(aerr, (size_t)k, sizeof(AtomError), cmp_atom_error_desc);

    /* how many to split */
    int budget = max_total_atoms - k;
    if (budget <= 0) { free(rendered); free(aerr); return 0; }
    int to_split = max_split;
    if (to_split > budget) to_split = budget;

    /* only split atoms that actually have error */
    int actual_split = 0;
    for (int s = 0; s < to_split; ++s) {
        if (aerr[s].error < 1e-6 || aerr[s].count < 2) break;
        actual_split++;
    }

    if (actual_split == 0) { free(rendered); free(aerr); return 0; }

    /* grow atom array */
    int new_k = k + actual_split;
    atoms = (ChootAtom*)realloc(atoms, (size_t)new_k * sizeof(ChootAtom));
    *p_atoms = atoms;

    for (int s = 0; s < actual_split; ++s) {
        int ci = aerr[s].idx;
        ChootAtom parent = atoms[ci];

        /* eigendecomposition of 2x2 covariance to find major axis */
        float sxx = parent.sxx, sxy = parent.sxy, syy = parent.syy;
        float trace = sxx + syy;
        float det   = sxx * syy - sxy * sxy;
        float disc  = trace * trace * 0.25f - det;
        if (disc < 0.0f) disc = 0.0f;
        float sq = sqrtf(disc);
        float lambda1 = trace * 0.5f + sq;  /* larger eigenvalue */

        /* eigenvector for lambda1 */
        float vx, vy;
        if (fabsf(sxy) > 1e-10f) {
            vx = lambda1 - syy;
            vy = sxy;
        } else {
            vx = (sxx >= syy) ? 1.0f : 0.0f;
            vy = (sxx >= syy) ? 0.0f : 1.0f;
        }
        float vlen = sqrtf(vx * vx + vy * vy);
        if (vlen > 1e-10f) { vx /= vlen; vy /= vlen; }

        /* offset along major axis by sqrt(eigenvalue) * 0.5 */
        float offset = sqrtf(lambda1) * 0.5f;

        /* child A (modify parent in place) */
        atoms[ci].x = clamp01(parent.x + vx * offset);
        atoms[ci].y = clamp01(parent.y + vy * offset);
        atoms[ci].sxx = parent.sxx * 0.5f;
        atoms[ci].sxy = parent.sxy * 0.5f;
        atoms[ci].syy = parent.syy * 0.5f;

        /* child B (new atom) */
        int ni = k + s;
        atoms[ni] = parent;
        atoms[ni].x = clamp01(parent.x - vx * offset);
        atoms[ni].y = clamp01(parent.y - vy * offset);
        atoms[ni].sxx = parent.sxx * 0.5f;
        atoms[ni].sxy = parent.sxy * 0.5f;
        atoms[ni].syy = parent.syy * 0.5f;
    }

    *p_k = new_k;

    free(rendered);
    free(aerr);
    return actual_split;
}

/* ------------------------------------------------------------------ */
/* Iterative analysis-by-synthesis colour correction                   */
/*                                                                     */
/* Unlike soft_refit (which fits each atom independently to the ground */
/* truth), this accounts for atom overlap:                             */
/*   1. Render all atoms.                                              */
/*   2. For each pixel, compute additive residual.                     */
/*   3. For each atom, compute Gaussian-weighted mean residual.        */
/*   4. Adjust atom colour by that residual (with damping).            */
/*                                                                     */
/* This converges because each step reduces the rendering error.       */
/* ------------------------------------------------------------------ */

static double abs_colour_correct(
    const Pixel5D* pixels,
    uint32_t img_w, uint32_t img_h,
    ChootAtom* atoms, int k,
    int num_iters, float damping)
{
    int n = (int)(img_w * img_h);
    float* rendered = (float*)malloc((size_t)n * 3 * sizeof(float));
    if (!rendered) return -1.0;

    /* Save best state so we can roll back if we diverge */
    ChootAtom* best_atoms = (ChootAtom*)malloc((size_t)k * sizeof(ChootAtom));
    if (!best_atoms) { free(rendered); return -1.0; }
    memcpy(best_atoms, atoms, (size_t)k * sizeof(ChootAtom));

    double best_mse = 1e30;

    int iter;
    for (iter = 0; iter < num_iters; ++iter) {
        /* render current state */
        ChootImage tmp;
        tmp.atom_count = (uint32_t)k;
        tmp.atoms = atoms;
        choot_render_linear_rgb(&tmp, (int)img_w, (int)img_h, rendered);

        /* compute MSE */
        double mse = 0.0;
        { int i;
        for (i = 0; i < n; ++i) {
            float dr = pixels[i].f[2] - rendered[i * 3 + 0];
            float dg = pixels[i].f[3] - rendered[i * 3 + 1];
            float db = pixels[i].f[4] - rendered[i * 3 + 2];
            mse += dr * dr + dg * dg + db * db;
        } }
        mse /= (n * 3);

        if (mse < best_mse) {
            best_mse = mse;
            memcpy(best_atoms, atoms, (size_t)k * sizeof(ChootAtom));
        }

        double psnr = (mse > 1e-10) ? 10.0 * log10(1.0 / mse) : 99.0;
        fprintf(stderr, "    iter %d: PSNR=%.2f dB\n", iter + 1, psnr);

        if (iter > 2 && mse > best_mse * 1.01) {
            fprintf(stderr, "    diverging — rolling back to best.\n");
            memcpy(atoms, best_atoms, (size_t)k * sizeof(ChootAtom));
            break;
        }

        /* for each atom, accumulate weighted residual and weighted GT */
        { int c;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16)
#endif
        for (c = 0; c < k; ++c) {
            if (atoms[c].alpha <= 0.0f) continue;
            ChootAtom* a = &atoms[c];
            float det = a->sxx * a->syy - a->sxy * a->sxy;
            if (det <= 0.0f) continue;
            float inv_det = 1.0f / det;
            float inv_xx =  a->syy * inv_det;
            float inv_xy = -a->sxy * inv_det;
            float inv_yy =  a->sxx * inv_det;

            /* bounding box: 4*sqrt(max eigenvalue) */
            float max_s = a->sxx > a->syy ? a->sxx : a->syy;
            float radius = 4.0f * sqrtf(max_s);
            int y0 = (int)((a->y - radius) * (float)img_h); if (y0 < 0) y0 = 0;
            int y1 = (int)((a->y + radius) * (float)img_h) + 1; if (y1 > (int)img_h) y1 = (int)img_h;
            int x0 = (int)((a->x - radius) * (float)img_w); if (x0 < 0) x0 = 0;
            int x1 = (int)((a->x + radius) * (float)img_w) + 1; if (x1 > (int)img_w) x1 = (int)img_w;

            /* accumulate weighted GT and weighted rendered */
            double wgt_r = 0, wgt_g = 0, wgt_b = 0;
            double wrd_r = 0, wrd_g = 0, wrd_b = 0;
            double w_sum = 0;

            int py, px;
            for (py = y0; py < y1; ++py) {
                for (px = x0; px < x1; ++px) {
                    int idx = py * (int)img_w + px;
                    float pxn = ((float)px + 0.5f) / (float)img_w;
                    float pyn = ((float)py + 0.5f) / (float)img_h;
                    float dx = pxn - a->x;
                    float dy = pyn - a->y;
                    float quad = dx * (inv_xx * dx + inv_xy * dy)
                               + dy * (inv_xy * dx + inv_yy * dy);
                    if (quad > 12.0f) continue;
                    float g = expf(-0.5f * quad);
                    float w = a->alpha * g;
                    if (w < 1e-7f) continue;

                    wgt_r += w * (double)pixels[idx].f[2];
                    wgt_g += w * (double)pixels[idx].f[3];
                    wgt_b += w * (double)pixels[idx].f[4];
                    wrd_r += w * (double)rendered[idx * 3 + 0];
                    wrd_g += w * (double)rendered[idx * 3 + 1];
                    wrd_b += w * (double)rendered[idx * 3 + 2];
                    w_sum += w;
                }
            }

            if (w_sum < 1e-8) continue;

            /* target = weighted average of GT under this atom's Gaussian */
            float tgt_r = (float)(wgt_r / w_sum);
            float tgt_g = (float)(wgt_g / w_sum);
            float tgt_b = (float)(wgt_b / w_sum);

            /* current rendered = weighted average of rendered under this atom */
            float cur_r = (float)(wrd_r / w_sum);
            float cur_g = (float)(wrd_g / w_sum);
            float cur_b = (float)(wrd_b / w_sum);

            /* residual: how much rendered differs from GT in this atom's region */
            float res_r = tgt_r - cur_r;
            float res_g = tgt_g - cur_g;
            float res_b = tgt_b - cur_b;

            /* get current atom color */
            float ar, ag, ab;
            choot_ycocg_to_rgb(a->Y, a->Co, a->Cg, &ar, &ag, &ab);

            /* adjust atom colour by damped residual */
            ar += damping * res_r;
            ag += damping * res_g;
            ab += damping * res_b;

            /* clamp */
            if (ar < 0) ar = 0; if (ar > 1) ar = 1;
            if (ag < 0) ag = 0; if (ag > 1) ag = 1;
            if (ab < 0) ab = 0; if (ab > 1) ab = 1;

            float luma, co, cg;
            choot_rgb_to_ycocg(ar, ag, ab, &luma, &co, &cg);
            a->Y  = luma;
            a->Co = co;
            a->Cg = cg;
        } }
    }

    /* ensure we end with the best state */
    memcpy(atoms, best_atoms, (size_t)k * sizeof(ChootAtom));
    free(best_atoms);
    free(rendered);
    return best_mse;
}

/* ------------------------------------------------------------------ */
/* Build 5-D centres from atom positions (for Lloyd re-clustering)     */
/* ------------------------------------------------------------------ */

static void atoms_to_centers(const ChootAtom* atoms, int k,
                             float spatial_weight, float* centers) {
    for (int c = 0; c < k; ++c) {
        centers[c * 5 + 0] = atoms[c].x * spatial_weight;
        centers[c * 5 + 1] = atoms[c].y * spatial_weight;
        float r, g, b;
        choot_ycocg_to_rgb(atoms[c].Y, atoms[c].Co, atoms[c].Cg, &r, &g, &b);
        centers[c * 5 + 2] = r;
        centers[c * 5 + 3] = g;
        centers[c * 5 + 4] = b;
    }
}

/* ------------------------------------------------------------------ */
/* Gradient-based color refinement (analysis-by-synthesis)              */
/*                                                                     */
/* For each optimization step:                                         */
/*   1. Render the current atom set at image resolution.               */
/*   2. Compute per-pixel residual = ground_truth - rendered.          */
/*   3. For each atom, accumulate the Gaussian-weighted residual.      */
/*   4. Adjust the atom's color in the direction that reduces error.   */
/*                                                                     */
/* This is much more effective than soft_refit because it accounts     */
/* for atom overlap: each atom adjusts based on the REMAINING error    */
/* after all other atoms contribute.                                   */
/* ------------------------------------------------------------------ */

static void gradient_color_refine(
    const Pixel5D* pixels, int n_pixels,
    uint32_t img_w, uint32_t img_h,
    ChootAtom* atoms, int k,
    int num_steps, float learning_rate)
{
    int n = (int)(img_w * img_h);
    float* rendered = (float*)malloc((size_t)n * 3 * sizeof(float));
    if (!rendered) return;

    int step;
    for (step = 0; step < num_steps; ++step) {
        /* render current state */
        ChootImage tmp;
        tmp.atom_count = (uint32_t)k;
        tmp.atoms = atoms;
        choot_render_linear_rgb(&tmp, (int)img_w, (int)img_h, rendered);

        /* for each atom, accumulate weighted residual */
        { int c;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16)
#endif
        for (c = 0; c < k; ++c) {
            if (atoms[c].alpha <= 0.0f) continue;
            ChootAtom* a = &atoms[c];
            float det = a->sxx * a->syy - a->sxy * a->sxy;
            if (det <= 0.0f) continue;
            float inv_det = 1.0f / det;
            float inv_xx =  a->syy * inv_det;
            float inv_xy = -a->sxy * inv_det;
            float inv_yy =  a->sxx * inv_det;

            double grad_r = 0, grad_g = 0, grad_b = 0, w_sum = 0;
            int py, px;

            /* bounding box from covariance: 4*sqrt(eigenvalue) radius */
            float max_r = 4.0f * sqrtf(a->sxx > a->syy ? a->sxx : a->syy);
            int y0 = (int)((a->y - max_r) * img_h); if (y0 < 0) y0 = 0;
            int y1 = (int)((a->y + max_r) * img_h) + 1; if (y1 > (int)img_h) y1 = (int)img_h;
            int x0 = (int)((a->x - max_r) * img_w); if (x0 < 0) x0 = 0;
            int x1 = (int)((a->x + max_r) * img_w) + 1; if (x1 > (int)img_w) x1 = (int)img_w;

            for (py = y0; py < y1; ++py) {
                for (px = x0; px < x1; ++px) {
                    int idx = py * (int)img_w + px;
                    float pxn = ((float)px + 0.5f) / (float)img_w;
                    float pyn = ((float)py + 0.5f) / (float)img_h;
                    float dx = pxn - a->x;
                    float dy = pyn - a->y;
                    float quad = dx * (inv_xx * dx + inv_xy * dy)
                               + dy * (inv_xy * dx + inv_yy * dy);
                    if (quad > 12.0f) continue;
                    float g = expf(-0.5f * quad);
                    float w = a->alpha * g;
                    if (w < 1e-6f) continue;

                    /* residual = ground_truth - rendered */
                    float dr = pixels[idx].f[2] - rendered[idx * 3 + 0];
                    float dg = pixels[idx].f[3] - rendered[idx * 3 + 1];
                    float db = pixels[idx].f[4] - rendered[idx * 3 + 2];

                    grad_r += w * dr;
                    grad_g += w * dg;
                    grad_b += w * db;
                    w_sum  += w;
                }
            }

            if (w_sum > 1e-8) {
                float lr = learning_rate;
                float cur_r, cur_g, cur_b;
                choot_ycocg_to_rgb(a->Y, a->Co, a->Cg, &cur_r, &cur_g, &cur_b);

                cur_r += lr * (float)(grad_r / w_sum);
                cur_g += lr * (float)(grad_g / w_sum);
                cur_b += lr * (float)(grad_b / w_sum);

                /* clamp to valid range */
                if (cur_r < 0) cur_r = 0; if (cur_r > 1) cur_r = 1;
                if (cur_g < 0) cur_g = 0; if (cur_g > 1) cur_g = 1;
                if (cur_b < 0) cur_b = 0; if (cur_b > 1) cur_b = 1;

                float luma, co, cg;
                choot_rgb_to_ycocg(cur_r, cur_g, cur_b, &luma, &co, &cg);
                a->Y  = luma;
                a->Co = co;
                a->Cg = cg;
            }
        } }
    }

    free(rendered);
}

/* ------------------------------------------------------------------ */
/* Recycle dead atoms into high-error regions                          */
/*                                                                     */
/* After re-clustering, empty clusters become dead (alpha=0).          */
/* Instead of throwing them away, move them to pixel locations with    */
/* the highest render error, giving them useful coverage.              */
/* ------------------------------------------------------------------ */

static void recycle_dead_atoms(
    const Pixel5D* pixels,
    uint32_t img_w, uint32_t img_h,
    ChootAtom* atoms, int k,
    float spatial_weight)
{
    int n = (int)(img_w * img_h);

    /* count dead atoms */
    int dead_count = 0;
    int ci;
    for (ci = 0; ci < k; ++ci)
        if (atoms[ci].alpha <= 0.0f) dead_count++;
    if (dead_count == 0) return;

    /* render and find highest-error pixels */
    float* rendered = (float*)malloc((size_t)n * 3 * sizeof(float));
    if (!rendered) return;
    ChootImage tmp;
    tmp.atom_count = (uint32_t)k;
    tmp.atoms = atoms;
    choot_render_linear_rgb(&tmp, (int)img_w, (int)img_h, rendered);

    /* compute per-pixel error */
    float* perr = (float*)malloc((size_t)n * sizeof(float));
    if (!perr) { free(rendered); return; }
    int pi;
    for (pi = 0; pi < n; ++pi) {
        float dr = pixels[pi].f[2] - rendered[pi * 3 + 0];
        float dg = pixels[pi].f[3] - rendered[pi * 3 + 1];
        float db = pixels[pi].f[4] - rendered[pi * 3 + 2];
        perr[pi] = dr*dr + dg*dg + db*db;
    }

    /* find top dead_count error pixels (simple selection) */
    int* top_idx = (int*)malloc((size_t)dead_count * sizeof(int));
    if (!top_idx) { free(perr); free(rendered); return; }
    int found = 0;

    /* partial sort: pick top dead_count by iterative max extraction */
    int max_to_find = dead_count;
    if (max_to_find > 500) max_to_find = 500; /* cap for perf */
    for (found = 0; found < max_to_find; ++found) {
        float best_err = -1.0f;
        int best_idx = 0;
        for (pi = 0; pi < n; ++pi) {
            if (perr[pi] > best_err) {
                best_err = perr[pi];
                best_idx = pi;
            }
        }
        if (best_err < 1e-6f) break;
        top_idx[found] = best_idx;
        /* suppress neighbours to spread atoms out */
        int bx = best_idx % (int)img_w;
        int by = best_idx / (int)img_w;
        int r = 3; /* suppression radius */
        int sy, sx;
        for (sy = by - r; sy <= by + r; ++sy) {
            if (sy < 0 || sy >= (int)img_h) continue;
            for (sx = bx - r; sx <= bx + r; ++sx) {
                if (sx < 0 || sx >= (int)img_w) continue;
                perr[sy * (int)img_w + sx] = 0.0f;
            }
        }
    }

    /* assign dead atoms to high-error positions */
    int placed = 0;
    float pix_sx = 1.5f / (float)img_w;
    float pix_sy = 1.5f / (float)img_h;
    for (ci = 0; ci < k && placed < found; ++ci) {
        if (atoms[ci].alpha > 0.0f) continue;
        int pidx = top_idx[placed];
        float px = ((float)(pidx % (int)img_w) + 0.5f) / (float)img_w;
        float py = ((float)(pidx / (int)img_w) + 0.5f) / (float)img_h;
        atoms[ci].x = px;
        atoms[ci].y = py;
        atoms[ci].sxx = pix_sx * pix_sx;
        atoms[ci].sxy = 0.0f;
        atoms[ci].syy = pix_sy * pix_sy;
        atoms[ci].alpha = 1.0f;

        float luma, co, cg;
        choot_rgb_to_ycocg(pixels[pidx].f[2], pixels[pidx].f[3], pixels[pidx].f[4],
                           &luma, &co, &cg);
        atoms[ci].Y  = luma;
        atoms[ci].Co = co;
        atoms[ci].Cg = cg;
        atoms[ci].flags = 0;
        placed++;
    }

    free(top_idx);
    free(perr);
    free(rendered);
    (void)spatial_weight;
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

static void usage(const char* exe) {
    fprintf(stderr,
        "Usage: %s input.png output.choot [options]\n"
        "\n"
        "Options:\n"
        "  atoms=N          Number of initial atoms         (default 5000)\n"
        "  max_atoms=N      Maximum atoms after splitting   (default 20000)\n"
        "  spatial=F        Position vs colour weight        (default 4.0)\n"
        "  sharp=F          Covariance sharpness (>1=sharper)(default 1.3)\n"
        "  refine=N         Number of refinement passes      (default 4)\n"
        "  split=N          Atoms to split per pass (pct)    (default 20)\n"
        "  iters=N          Lloyd iterations per stage       (default 40)\n"
        "  sample=N         Pixel sampling stride (1=full)   (default 2)\n"
        "  fast             Speed preset (lower quality)\n"
        "  hq               Quality preset (slower)\n"
        "\n"
        "Examples:\n"
        "  %s photo.png photo.choot                        (defaults: good quality)\n"
        "  %s photo.png hq.choot hq                        (highest quality)\n"
        "  %s photo.png fast.choot fast                   (fast preview)\n"
        "  %s icon.png  icon.choot atoms=500               (small file)\n",
        exe, exe, exe, exe, exe);
}

static int parse_int_opt(const char* arg, const char* name) {
    size_t nlen = strlen(name);
    if (strncmp(arg, name, nlen) == 0 && arg[nlen] == '=')
        return atoi(arg + nlen + 1);
    return -1;
}

static float parse_float_opt(const char* arg, const char* name) {
    size_t nlen = strlen(name);
    if (strncmp(arg, name, nlen) == 0 && arg[nlen] == '=')
        return (float)atof(arg + nlen + 1);
    return -1.0f;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        usage(argv[0]);
        return 1;
    }

    const char* in_path  = argv[1];
    const char* out_path = argv[2];

    /* defaults tuned for quality */
    int   init_atoms     = 5000;
    int   max_atoms      = 40000;
    float spatial_weight = 2.0f;
    float sharpness      = 1.2f;
    int   refine_passes  = 6;
    int   split_pct      = 30;
    int   lloyd_iters    = 50;
    int   sample_stride  = 2;
    int   grad_steps     = 0;
    int   fast_preset    = 0;
    int   hq_preset      = 0;

    for (int i = 3; i < argc; ++i) {
        int v;
        float fv;
        if ((v = parse_int_opt(argv[i], "atoms"))      > 0) init_atoms     = v;
        if ((v = parse_int_opt(argv[i], "max_atoms"))   > 0) max_atoms      = v;
        if ((fv = parse_float_opt(argv[i], "spatial"))  > 0) spatial_weight  = fv;
        if ((fv = parse_float_opt(argv[i], "sharp"))    > 0) sharpness      = fv;
        if ((v = parse_int_opt(argv[i], "refine"))     >= 0) refine_passes  = v;
        if ((v = parse_int_opt(argv[i], "split"))       > 0) split_pct      = v;
        if ((v = parse_int_opt(argv[i], "iters"))       > 0) lloyd_iters    = v;
        if ((v = parse_int_opt(argv[i], "sample"))      > 0) sample_stride  = v;
        if (strcmp(argv[i], "fast") == 0) fast_preset = 1;
        if (strcmp(argv[i], "hq") == 0)   hq_preset   = 1;
    }

    if (fast_preset) {
        init_atoms     = 2000;
        max_atoms      = 8000;
        spatial_weight = 2.0f;
        sharpness      = 1.1f;
        refine_passes  = 2;
        split_pct      = 20;
        lloyd_iters    = 20;
        sample_stride  = 3;
        grad_steps     = 0;
    }
    if (hq_preset) {
        init_atoms     = 8000;
        max_atoms      = 46000;
        spatial_weight = 2.0f;
        sharpness      = 1.3f;
        refine_passes  = 8;
        split_pct      = 30;
        lloyd_iters    = 60;
        sample_stride  = 1;
        grad_steps     = 0;
    }

    if (max_atoms < init_atoms) max_atoms = init_atoms;

    fprintf(stderr, "CHOOT encoder v2 (multi-stage)\n");
    fprintf(stderr, "  init_atoms = %d\n", init_atoms);
    fprintf(stderr, "  max_atoms  = %d\n", max_atoms);
    fprintf(stderr, "  spatial    = %.1f\n", spatial_weight);
    fprintf(stderr, "  sharp      = %.2f\n", sharpness);
    fprintf(stderr, "  refine     = %d passes\n", refine_passes);
    fprintf(stderr, "  split      = %d%% per pass\n", split_pct);
    fprintf(stderr, "  iters      = %d per stage\n", lloyd_iters);
    fprintf(stderr, "  sample     = %d\n", sample_stride);

    /* ---- load PNG ---- */
    WicImage img;
    char err[256];
    if (!wic_load_png(in_path, &img, err, sizeof(err))) {
        fprintf(stderr, "PNG load failed: %s\n", err);
        return 1;
    }
    fprintf(stderr, "Loaded %ux%u PNG (%u pixels)\n",
            img.width, img.height, img.width * img.height);

    int n = (int)(img.width * img.height);
    if (init_atoms > n) init_atoms = n;

    /* ---- build 5-D features ---- */
    Pixel5D* pixels = (Pixel5D*)malloc((size_t)n * sizeof(Pixel5D));
    if (!pixels) { fprintf(stderr, "OOM\n"); wic_free(&img); return 1; }

    for (int py = 0; py < (int)img.height; ++py) {
        for (int px = 0; px < (int)img.width; ++px) {
            int idx = py * (int)img.width + px;
            const uint8_t* src = img.rgba + (size_t)idx * 4;
            float nx = ((float)px + 0.5f) / (float)img.width;
            float ny = ((float)py + 0.5f) / (float)img.height;
            pixels[idx].f[0] = nx * spatial_weight;
            pixels[idx].f[1] = ny * spatial_weight;
            pixels[idx].f[2] = srgb_to_linear(src[0] / 255.0f);
            pixels[idx].f[3] = srgb_to_linear(src[1] / 255.0f);
            pixels[idx].f[4] = srgb_to_linear(src[2] / 255.0f);
            pixels[idx].a    = src[3] / 255.0f;
        }
    }

    /* ==== Stage 1: initial k-means ==== */
    fprintf(stderr, "\n[Stage 1] k-means++ seeding %d atoms...\n", init_atoms);
    int k = init_atoms;
    float* centers  = (float*)calloc((size_t)k * 5, sizeof(float));
    float* dist_buf = (float*)malloc((size_t)n * sizeof(float));
    int*   assign   = (int*)calloc((size_t)n, sizeof(int));

    kmeans_pp_seed(pixels, n, centers, k, dist_buf);
    free(dist_buf);

    fprintf(stderr, "[Stage 1] Lloyd refining (%d iters)...\n", lloyd_iters);
    kmeans_lloyd(pixels, n, centers, k, assign, lloyd_iters);

    fprintf(stderr, "[Stage 1] Fitting initial atoms...\n");
    ChootAtom* atoms = (ChootAtom*)calloc((size_t)k, sizeof(ChootAtom));
    fit_atoms(pixels, n, assign, k, spatial_weight, sharpness, img.width, img.height, atoms);

    fprintf(stderr, "[Stage 1] Soft colour refit...\n");
    soft_refit_colours(pixels, n, img.width, img.height, atoms, k, spatial_weight, sample_stride);

    /* ==== Stage 2: iterative split + refit ==== */
    for (int pass = 0; pass < refine_passes; ++pass) {
        int to_split = (k * split_pct) / 100;
        if (to_split < 10) to_split = 10;
        if (to_split > 3000) to_split = 3000;

        fprintf(stderr, "\n[Stage 2, pass %d/%d] Splitting up to %d atoms (current: %d)...\n",
                pass + 1, refine_passes, to_split, k);

        int did_split = split_worst_atoms(
            pixels, n, img.width, img.height,
            &atoms, &k, to_split, max_atoms, spatial_weight, sample_stride);

        if (did_split == 0) {
            fprintf(stderr, "  No more atoms to split. Done.\n");
            break;
        }
        fprintf(stderr, "  Split %d -> now %d total\n", did_split, k);

        fprintf(stderr, "  Soft colour refit...\n");
        soft_refit_colours(pixels, n, img.width, img.height, atoms, k, spatial_weight, sample_stride);
    }

    /* ==== Stage 3: final polish ==== */
    fprintf(stderr, "\n[Stage 3] Final soft colour refit (3 passes at full res)...\n");
    soft_refit_colours(pixels, n, img.width, img.height, atoms, k, spatial_weight, 1);
    soft_refit_colours(pixels, n, img.width, img.height, atoms, k, spatial_weight, 1);
    soft_refit_colours(pixels, n, img.width, img.height, atoms, k, spatial_weight, 1);

    /* ---- compute final quality metric ---- */
    {
        float* rendered = (float*)malloc((size_t)img.width * img.height * 3 * sizeof(float));
        if (rendered) {
            ChootImage tmp;
            tmp.atom_count = (uint32_t)k;
            tmp.atoms = atoms;
            choot_render_linear_rgb(&tmp, (int)img.width, (int)img.height, rendered);

            double mse = 0.0;
            for (int i = 0; i < n; ++i) {
                float dr = pixels[i].f[2] - rendered[i * 3 + 0];
                float dg = pixels[i].f[3] - rendered[i * 3 + 1];
                float db = pixels[i].f[4] - rendered[i * 3 + 2];
                mse += dr * dr + dg * dg + db * db;
            }
            mse /= (n * 3);
            double psnr = (mse > 1e-10) ? 10.0 * log10(1.0 / mse) : 99.0;
            fprintf(stderr, "\n  Quality: MSE=%.6f  PSNR=%.2f dB\n", mse, psnr);
            free(rendered);
        }
    }

    free(centers);
    free(assign);
    free(pixels);
    wic_free(&img);

    /* ---- prune dead atoms ---- */
    int alive = 0;
    for (int i = 0; i < k; ++i)
        if (atoms[i].alpha > 0.0f) alive++;

    if (alive < k) {
        ChootAtom* compacted = (ChootAtom*)malloc((size_t)alive * sizeof(ChootAtom));
        int j = 0;
        for (int i = 0; i < k; ++i)
            if (atoms[i].alpha > 0.0f) compacted[j++] = atoms[i];
        free(atoms);
        atoms = compacted;
        k = alive;
    }

    /* ---- write .choot ---- */
    ChootImage out;
    out.atom_count = (uint32_t)k;
    out.atoms = atoms;

    if (!choot_save(out_path, &out, err, sizeof(err))) {
        fprintf(stderr, "Save failed: %s\n", err);
        free(atoms);
        return 1;
    }

    fprintf(stderr, "\nWrote %d atoms -> %s\n", k, out_path);
    fprintf(stderr, "File size: %u bytes (%.1f KB)\n",
            (unsigned)(CHOOT_HEADER_SIZE + k * 20),
            (float)(CHOOT_HEADER_SIZE + k * 20) / 1024.0f);

    free(atoms);
    return 0;
}
