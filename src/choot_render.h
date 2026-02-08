#pragma once

#include <stdint.h>
#include "choot.h"

#ifdef __cplusplus
extern "C" {
#endif

int choot_render_linear_rgb(const ChootImage* image, int width, int height, float* out_rgb);
int choot_render_rgba8(const ChootImage* image, int width, int height, uint8_t* out_rgba);

#ifdef __cplusplus
}
#endif
