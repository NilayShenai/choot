#pragma once

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CHOOT_MAGIC "CHOOT\0\0"
#define CHOOT_HEADER_SIZE 24
#define CHOOT_VERSION 0

typedef struct ChootAtom {
    float x;
    float y;
    float sxx;
    float sxy;
    float syy;
    float alpha;
    float Y;
    float Co;
    float Cg;
    uint16_t flags;
} ChootAtom;

typedef struct ChootImage {
    uint32_t atom_count;
    ChootAtom* atoms;
} ChootImage;

int choot_load(const char* path, ChootImage* out_image, char* err, size_t err_cap);
int choot_save(const char* path, const ChootImage* image, char* err, size_t err_cap);
void choot_free(ChootImage* image);

void choot_rgb_to_ycocg(float r, float g, float b, float* out_y, float* out_co, float* out_cg);
void choot_ycocg_to_rgb(float y, float co, float cg, float* out_r, float* out_g, float* out_b);

#ifdef __cplusplus
}
#endif
