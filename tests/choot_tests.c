#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdint.h>

#include "choot.h"
#include "choot_render.h"

static int g_failures = 0;

#define ASSERT_TRUE(expr, msg) do { \
    if (!(expr)) { \
        fprintf(stderr, "[FAIL] %s (%s:%d)\n", msg, __FILE__, __LINE__); \
        g_failures++; \
    } \
} while (0)

#define ASSERT_NEAR(a, b, eps, msg) do { \
    float _da = (float)(a); \
    float _db = (float)(b); \
    float _diff = fabsf(_da - _db); \
    if (_diff > (eps)) { \
        fprintf(stderr, "[FAIL] %s (%.6f vs %.6f, eps=%.6f)\n", msg, _da, _db, (float)(eps)); \
        g_failures++; \
    } \
} while (0)

static int write_file(const char* path, const void* data, size_t size) {
    FILE* f = fopen(path, "wb");
    if (!f) return 0;
    size_t w = fwrite(data, 1, size, f);
    fclose(f);
    return w == size;
}

static int file_exists(const char* path) {
    FILE* f = fopen(path, "rb");
    if (!f) return 0;
    fclose(f);
    return 1;
}

static void make_red_atom(ChootAtom* a, float x, float y, float sigma) {
    float Y, Co, Cg;
    choot_rgb_to_ycocg(1.0f, 0.0f, 0.0f, &Y, &Co, &Cg);
    a->x = x;
    a->y = y;
    a->sxx = sigma * sigma;
    a->sxy = 0.0f;
    a->syy = sigma * sigma;
    a->alpha = 1.0f;
    a->Y = Y;
    a->Co = Co;
    a->Cg = Cg;
    a->flags = 0;
}

static void test_invalid_magic(void) {
    unsigned char header[CHOOT_HEADER_SIZE];
    memset(header, 0, sizeof(header));
    memcpy(header, "BADMAGC", 7);
    header[8] = 0;
    header[9] = 0;
    header[16] = CHOOT_HEADER_SIZE;

    const char* path = "test_invalid_magic.choot";
    ASSERT_TRUE(write_file(path, header, sizeof(header)), "write invalid magic");

    ChootImage img;
    char err[128] = {0};
    int ok = choot_load(path, &img, err, sizeof(err));
    ASSERT_TRUE(!ok, "invalid magic rejected");
    ASSERT_TRUE(err[0] != 0, "error string set");
    remove(path);
}

static void test_invalid_version(void) {
    unsigned char header[CHOOT_HEADER_SIZE];
    memset(header, 0, sizeof(header));
    memcpy(header, CHOOT_MAGIC, 8);
    header[8] = 1; /* version = 1 */
    header[9] = 0;
    header[16] = CHOOT_HEADER_SIZE;

    const char* path = "test_invalid_version.choot";
    ASSERT_TRUE(write_file(path, header, sizeof(header)), "write invalid version");

    ChootImage img;
    char err[128] = {0};
    int ok = choot_load(path, &img, err, sizeof(err));
    ASSERT_TRUE(!ok, "invalid version rejected");
    ASSERT_TRUE(err[0] != 0, "error string set");
    remove(path);
}

static void test_invalid_header_size(void) {
    unsigned char header[CHOOT_HEADER_SIZE];
    memset(header, 0, sizeof(header));
    memcpy(header, CHOOT_MAGIC, 8);
    header[8] = 0;
    header[9] = 0;
    header[16] = 32; /* header_size != 24 */

    const char* path = "test_invalid_header.choot";
    ASSERT_TRUE(write_file(path, header, sizeof(header)), "write invalid header size");

    ChootImage img;
    char err[128] = {0};
    int ok = choot_load(path, &img, err, sizeof(err));
    ASSERT_TRUE(!ok, "invalid header size rejected");
    ASSERT_TRUE(err[0] != 0, "error string set");
    remove(path);
}

static void test_truncated_file(void) {
    unsigned char header[CHOOT_HEADER_SIZE];
    memset(header, 0, sizeof(header));
    memcpy(header, CHOOT_MAGIC, 8);
    header[8] = 0;
    header[9] = 0;
    header[12] = 1; /* atom_count = 1 */
    header[16] = CHOOT_HEADER_SIZE;

    const char* path = "test_truncated.choot";
    ASSERT_TRUE(write_file(path, header, sizeof(header)), "write truncated file");

    ChootImage img;
    char err[128] = {0};
    int ok = choot_load(path, &img, err, sizeof(err));
    ASSERT_TRUE(!ok, "truncated file rejected");
    ASSERT_TRUE(err[0] != 0, "error string set");
    remove(path);
}

static void test_determinism(void) {
    ChootAtom atom;
    make_red_atom(&atom, 0.5f, 0.5f, 0.15f);
    ChootImage img;
    img.atom_count = 1;
    img.atoms = &atom;

    int w = 32, h = 32;
    size_t count = (size_t)w * (size_t)h * 3;
    float* a = (float*)malloc(count * sizeof(float));
    float* b = (float*)malloc(count * sizeof(float));

    ASSERT_TRUE(choot_render_linear_rgb(&img, w, h, a), "render A");
    ASSERT_TRUE(choot_render_linear_rgb(&img, w, h, b), "render B");

    int identical = memcmp(a, b, count * sizeof(float)) == 0;
    ASSERT_TRUE(identical, "deterministic render");

    free(a);
    free(b);
}

static void test_visual_symmetry(void) {
    ChootAtom atom;
    make_red_atom(&atom, 0.5f, 0.5f, 0.12f);
    ChootImage img;
    img.atom_count = 1;
    img.atoms = &atom;

    int w = 5, h = 5;
    float buf[5 * 5 * 3];
    ASSERT_TRUE(choot_render_linear_rgb(&img, w, h, buf), "render symmetry");

    float center = buf[(2 * w + 2) * 3 + 0];
    float edge   = buf[(0 * w + 2) * 3 + 0];
    float corner = buf[(0 * w + 0) * 3 + 0];

    ASSERT_TRUE(center > edge, "center > edge");
    ASSERT_TRUE(edge > corner, "edge > corner");
}

static void test_color_correctness(void) {
    ChootAtom atom;
    make_red_atom(&atom, 0.5f, 0.5f, 0.25f);
    ChootImage img;
    img.atom_count = 1;
    img.atoms = &atom;

    float buf[3];
    ASSERT_TRUE(choot_render_linear_rgb(&img, 1, 1, buf), "render 1x1");

    ASSERT_TRUE(buf[0] > 0.9f, "red channel high");
    ASSERT_TRUE(buf[1] < 0.1f, "green channel low");
    ASSERT_TRUE(buf[2] < 0.1f, "blue channel low");
}

static void test_stress_many_atoms(void) {
    const int k = 2000;
    ChootAtom* atoms = (ChootAtom*)calloc((size_t)k, sizeof(ChootAtom));
    ASSERT_TRUE(atoms != NULL, "alloc atoms");

    for (int i = 0; i < k; ++i) {
        float fx = (float)(i % 100) / 99.0f;
        float fy = (float)(i / 100) / 19.0f;
        atoms[i].x = fx;
        atoms[i].y = fy;
        atoms[i].sxx = 0.0003f;
        atoms[i].sxy = 0.0f;
        atoms[i].syy = 0.0003f;
        atoms[i].alpha = 0.25f;
        atoms[i].Y = 0.5f;
        atoms[i].Co = 0.0f;
        atoms[i].Cg = 0.0f;
        atoms[i].flags = 0;
    }

    ChootImage img;
    img.atom_count = k;
    img.atoms = atoms;

    int w = 64, h = 64;
    size_t count = (size_t)w * (size_t)h * 3;
    float* buf = (float*)malloc(count * sizeof(float));
    ASSERT_TRUE(choot_render_linear_rgb(&img, w, h, buf), "render many atoms");

    for (size_t i = 0; i < count; ++i) {
        ASSERT_TRUE(isfinite(buf[i]), "finite output");
    }

    free(buf);
    free(atoms);
}

static void test_tiny_atoms(void) {
    ChootAtom atom;
    make_red_atom(&atom, 0.5f, 0.5f, 1e-4f);
    ChootImage img;
    img.atom_count = 1;
    img.atoms = &atom;

    int w = 32, h = 32;
    size_t count = (size_t)w * (size_t)h * 3;
    float* buf = (float*)malloc(count * sizeof(float));
    ASSERT_TRUE(choot_render_linear_rgb(&img, w, h, buf), "render tiny atom");

    for (size_t i = 0; i < count; ++i) {
        ASSERT_TRUE(isfinite(buf[i]), "finite output tiny");
    }

    free(buf);
}

static void test_save_load_roundtrip(void) {
    ChootAtom atom;
    make_red_atom(&atom, 0.33f, 0.66f, 0.2f);
    ChootImage img;
    img.atom_count = 1;
    img.atoms = &atom;

    const char* path = "test_roundtrip.choot";
    char err[128] = {0};
    ASSERT_TRUE(choot_save(path, &img, err, sizeof(err)), "save roundtrip");
    ASSERT_TRUE(file_exists(path), "file exists");

    ChootImage loaded;
    ASSERT_TRUE(choot_load(path, &loaded, err, sizeof(err)), "load roundtrip");
    ASSERT_TRUE(loaded.atom_count == 1, "roundtrip atom_count");

    choot_free(&loaded);
    remove(path);
}

int main(void) {
    fprintf(stderr, "CHOOT tests\n");

    test_invalid_magic();
    test_invalid_version();
    test_invalid_header_size();
    test_truncated_file();

    test_determinism();
    test_visual_symmetry();
    test_color_correctness();

    test_stress_many_atoms();
    test_tiny_atoms();
    test_save_load_roundtrip();

    if (g_failures == 0) {
        fprintf(stderr, "All tests passed.\n");
        return 0;
    }

    fprintf(stderr, "%d tests failed.\n", g_failures);
    return 1;
}
