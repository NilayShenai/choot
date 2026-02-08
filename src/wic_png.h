#pragma once

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct WicImage {
    uint32_t width;
    uint32_t height;
    uint8_t* rgba; 
} WicImage;

int wic_load_png(const char* path, WicImage* out_image, char* err, size_t err_cap);
void wic_free(WicImage* image);

#ifdef __cplusplus
}
#endif
