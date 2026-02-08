#include "choot_render.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>

static float clamp01(float v) {
    if (v < 0.0f) return 0.0f;
    if (v > 1.0f) return 1.0f;
    return v;
}

static float linear_to_srgb(float c) {
    if (c <= 0.0031308f) return 12.92f * c;
    return 1.055f * powf(c, 1.0f / 2.4f) - 0.055f;
}

int choot_render_linear_rgb(const ChootImage* image, int width, int height, float* out_rgb) {
    if (!image || !out_rgb || width <= 0 || height <= 0) {
        return 0;
    }

    memset(out_rgb, 0, (size_t)width * (size_t)height * 3 * sizeof(float));

    { int y;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            float px = (x + 0.5f) / (float)width;
            float py = (y + 0.5f) / (float)height;

            float acc_r = 0.0f;
            float acc_g = 0.0f;
            float acc_b = 0.0f;
            float acc_a = 0.0f;

            for (uint32_t i = 0; i < image->atom_count; ++i) {
                const ChootAtom* a = &image->atoms[i];
                float dx = px - a->x;
                float dy = py - a->y;

                float det = a->sxx * a->syy - a->sxy * a->sxy;
                if (det <= 0.0f) {
                    continue;
                }
                float inv_det = 1.0f / det;
                float inv_xx = a->syy * inv_det;
                float inv_xy = -a->sxy * inv_det;
                float inv_yy = a->sxx * inv_det;

                float quad = dx * (inv_xx * dx + inv_xy * dy) + dy * (inv_xy * dx + inv_yy * dy);
                float g = expf(-0.5f * quad);
                float w = a->alpha * g;
                if (w <= 0.0f) {
                    continue;
                }

                float r, gcol, b;
                choot_ycocg_to_rgb(a->Y, a->Co, a->Cg, &r, &gcol, &b);

                acc_r += w * r;
                acc_g += w * gcol;
                acc_b += w * b;
                acc_a += w;
            }

            float out_r = 0.0f;
            float out_g = 0.0f;
            float out_b = 0.0f;
            if (acc_a > 1e-8f) {
                out_r = acc_r / acc_a;
                out_g = acc_g / acc_a;
                out_b = acc_b / acc_a;
            }

            size_t idx = ((size_t)y * (size_t)width + (size_t)x) * 3;
            out_rgb[idx + 0] = out_r;
            out_rgb[idx + 1] = out_g;
            out_rgb[idx + 2] = out_b;
        }
    } }

    return 1;
}

int choot_render_rgba8(const ChootImage* image, int width, int height, uint8_t* out_rgba) {
    if (!image || !out_rgba || width <= 0 || height <= 0) {
        return 0;
    }

    size_t pixel_count = (size_t)width * (size_t)height;
    float* linear = (float*)malloc(pixel_count * 3 * sizeof(float));
    if (!linear) {
        return 0;
    }

    if (!choot_render_linear_rgb(image, width, height, linear)) {
        free(linear);
        return 0;
    }

    { ptrdiff_t i;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (i = 0; i < (ptrdiff_t)pixel_count; ++i) {
        float out_r = clamp01(linear[i * 3 + 0]);
        float out_g = clamp01(linear[i * 3 + 1]);
        float out_b = clamp01(linear[i * 3 + 2]);

        float sr = linear_to_srgb(out_r);
        float sg = linear_to_srgb(out_g);
        float sb = linear_to_srgb(out_b);

        out_rgba[i * 4 + 0] = (uint8_t)(clamp01(sr) * 255.0f + 0.5f);
        out_rgba[i * 4 + 1] = (uint8_t)(clamp01(sg) * 255.0f + 0.5f);
        out_rgba[i * 4 + 2] = (uint8_t)(clamp01(sb) * 255.0f + 0.5f);
        out_rgba[i * 4 + 3] = 255;
    } }

    free(linear);
    return 1;
}
