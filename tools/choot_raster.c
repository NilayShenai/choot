#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "choot.h"
#include "choot_render.h"

static int write_ppm(const char* path, int width, int height, const uint8_t* rgba, int ascii_p3) {
    FILE* f = fopen(path, "wb");
    if (!f) return 0;

    if (!ascii_p3) {
        fprintf(f, "P6\n%d %d\n255\n", width, height);
        for (int y = 0; y < height; ++y) {
            const uint8_t* row = rgba + (size_t)y * (size_t)width * 4;
            for (int x = 0; x < width; ++x) {
                fwrite(row + x * 4, 1, 3, f);
            }
        }
    } else {
        fprintf(f, "P3\n%d %d\n255\n", width, height);
        for (int y = 0; y < height; ++y) {
            const uint8_t* row = rgba + (size_t)y * (size_t)width * 4;
            for (int x = 0; x < width; ++x) {
                const uint8_t* px = row + x * 4;
                fprintf(f, "%u %u %u\n", (unsigned)px[0], (unsigned)px[1], (unsigned)px[2]);
            }
        }
    }

    fclose(f);
    return 1;
}

int main(int argc, char** argv) {
    if (argc < 5) {
        fprintf(stderr, "Usage: %s input.choot output.ppm width height [--p3]\n", argv[0]);
        return 1;
    }

    const char* in_path = argv[1];
    const char* out_path = argv[2];
    int width = atoi(argv[3]);
    int height = atoi(argv[4]);

    if (width <= 0 || height <= 0) {
        fprintf(stderr, "Invalid resolution\n");
        return 1;
    }

    ChootImage img;
    char err[256];
    if (!choot_load(in_path, &img, err, sizeof(err))) {
        fprintf(stderr, "Load failed: %s\n", err);
        return 1;
    }

    size_t size = (size_t)width * (size_t)height * 4;
    uint8_t* rgba = (uint8_t*)malloc(size);
    if (!rgba) {
        fprintf(stderr, "Out of memory\n");
        choot_free(&img);
        return 1;
    }

    if (!choot_render_rgba8(&img, width, height, rgba)) {
        fprintf(stderr, "Render failed\n");
        free(rgba);
        choot_free(&img);
        return 1;
    }

    int ascii_p3 = 0;
    if (argc >= 6 && strcmp(argv[5], "--p3") == 0) {
        ascii_p3 = 1;
    }

    if (!write_ppm(out_path, width, height, rgba, ascii_p3)) {
        fprintf(stderr, "Write failed\n");
        free(rgba);
        choot_free(&img);
        return 1;
    }

    free(rgba);
    choot_free(&img);
    return 0;
}
