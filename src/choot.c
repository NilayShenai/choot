#include "choot.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void set_err(char* err, size_t cap, const char* msg) {
    if (!err || cap == 0) {
        return;
    }
    strncpy(err, msg, cap - 1);
    err[cap - 1] = '\0';
}

static uint16_t f32_to_f16(float f) {
    union { float f; uint32_t u; } in;
    in.f = f;
    uint32_t sign = (in.u >> 16) & 0x8000u;
    uint32_t mantissa = in.u & 0x007FFFFFu;
    uint32_t exp = in.u & 0x7F800000u;

    if (exp >= 0x47800000u) { 
        if (exp == 0x7F800000u && mantissa != 0) {
            return (uint16_t)(sign | 0x7E00u); 
        }
        return (uint16_t)(sign | 0x7C00u); 
    }

    if (exp <= 0x38000000u) {
        if (exp < 0x33000000u) {
            return (uint16_t)sign; 
        }
        uint32_t shifted = (mantissa | 0x00800000u) >> (113 - (exp >> 23));
        return (uint16_t)(sign | (shifted + 0x00001000u) >> 13);
    }

    uint16_t out = (uint16_t)(sign | ((exp - 0x38000000u) >> 13) | ((mantissa + 0x00001000u) >> 13));
    return out;
}

static float f16_to_f32(uint16_t h) {
    uint32_t sign = (uint32_t)(h & 0x8000u) << 16;
    uint32_t exp = (h >> 10) & 0x1Fu;
    uint32_t mantissa = h & 0x03FFu;
    uint32_t out;

    if (exp == 0) {
        if (mantissa == 0) {
            out = sign;
        } else {
            exp = 1;
            while ((mantissa & 0x0400u) == 0) {
                mantissa <<= 1;
                exp--;
            }
            mantissa &= 0x03FFu;
            out = sign | ((exp + 127 - 15) << 23) | (mantissa << 13);
        }
    } else if (exp == 31) {
        out = sign | 0x7F800000u | (mantissa << 13);
    } else {
        out = sign | ((exp + 127 - 15) << 23) | (mantissa << 13);
    }

    union { uint32_t u; float f; } conv;
    conv.u = out;
    return conv.f;
}

void choot_rgb_to_ycocg(float r, float g, float b, float* out_y, float* out_co, float* out_cg) {
    /* Standard YCoCg: lossless round-trip with the inverse below */
    float y  =  0.25f * r + 0.5f * g + 0.25f * b;
    float co =  0.5f  * r            - 0.5f  * b;
    float cg = -0.25f * r + 0.5f * g - 0.25f * b;
    if (out_y) *out_y = y;
    if (out_co) *out_co = co;
    if (out_cg) *out_cg = cg;
}

void choot_ycocg_to_rgb(float y, float co, float cg, float* out_r, float* out_g, float* out_b) {
    /* Exact inverse of choot_rgb_to_ycocg */
    float r = y + co - cg;
    float g = y + cg;
    float b = y - co - cg;
    if (out_r) *out_r = r;
    if (out_g) *out_g = g;
    if (out_b) *out_b = b;
}

int choot_load(const char* path, ChootImage* out_image, char* err, size_t err_cap) {
    if (!out_image) {
        set_err(err, err_cap, "null out_image");
        return 0;
    }
    memset(out_image, 0, sizeof(*out_image));

    FILE* f = fopen(path, "rb");
    if (!f) {
        set_err(err, err_cap, "failed to open file");
        return 0;
    }

    unsigned char header[CHOOT_HEADER_SIZE];
    if (fread(header, 1, CHOOT_HEADER_SIZE, f) != CHOOT_HEADER_SIZE) {
        fclose(f);
        set_err(err, err_cap, "failed to read header");
        return 0;
    }

    if (memcmp(header, CHOOT_MAGIC, 8) != 0) {
        fclose(f);
        set_err(err, err_cap, "bad magic");
        return 0;
    }

    uint16_t version = (uint16_t)(header[8] | (header[9] << 8));
    uint16_t flags = (uint16_t)(header[10] | (header[11] << 8));
    uint32_t atom_count = (uint32_t)(header[12] | (header[13] << 8) | (header[14] << 16) | (header[15] << 24));
    uint32_t header_size = (uint32_t)(header[16] | (header[17] << 8) | (header[18] << 16) | (header[19] << 24));

    if (version != CHOOT_VERSION) {
        fclose(f);
        set_err(err, err_cap, "unsupported version");
        return 0;
    }
    if (flags != 0) {
        fclose(f);
        set_err(err, err_cap, "unsupported flags");
        return 0;
    }
    if (header_size != CHOOT_HEADER_SIZE) {
        fclose(f);
        set_err(err, err_cap, "unsupported header size");
        return 0;
    }

    if (atom_count > 10000000u) {
        fclose(f);
        set_err(err, err_cap, "atom_count too large");
        return 0;
    }

    ChootAtom* atoms = (ChootAtom*)calloc(atom_count, sizeof(ChootAtom));
    if (!atoms) {
        fclose(f);
        set_err(err, err_cap, "out of memory");
        return 0;
    }

    for (uint32_t i = 0; i < atom_count; ++i) {
        unsigned char rec[20];
        if (fread(rec, 1, sizeof(rec), f) != sizeof(rec)) {
            free(atoms);
            fclose(f);
            set_err(err, err_cap, "failed to read atom");
            return 0;
        }

        uint16_t x = (uint16_t)(rec[0] | (rec[1] << 8));
        uint16_t y = (uint16_t)(rec[2] | (rec[3] << 8));
        uint16_t sxx = (uint16_t)(rec[4] | (rec[5] << 8));
        uint16_t sxy = (uint16_t)(rec[6] | (rec[7] << 8));
        uint16_t syy = (uint16_t)(rec[8] | (rec[9] << 8));
        uint16_t alpha = (uint16_t)(rec[10] | (rec[11] << 8));
        uint16_t Y = (uint16_t)(rec[12] | (rec[13] << 8));
        uint16_t Co = (uint16_t)(rec[14] | (rec[15] << 8));
        uint16_t Cg = (uint16_t)(rec[16] | (rec[17] << 8));
        uint16_t flags_atom = (uint16_t)(rec[18] | (rec[19] << 8));

        atoms[i].x = f16_to_f32(x);
        atoms[i].y = f16_to_f32(y);
        atoms[i].sxx = f16_to_f32(sxx);
        atoms[i].sxy = f16_to_f32(sxy);
        atoms[i].syy = f16_to_f32(syy);
        atoms[i].alpha = f16_to_f32(alpha);
        atoms[i].Y = f16_to_f32(Y);
        atoms[i].Co = f16_to_f32(Co);
        atoms[i].Cg = f16_to_f32(Cg);
        atoms[i].flags = flags_atom;
    }

    fclose(f);
    out_image->atom_count = atom_count;
    out_image->atoms = atoms;
    return 1;
}

int choot_save(const char* path, const ChootImage* image, char* err, size_t err_cap) {
    if (!image || (image->atom_count > 0 && !image->atoms)) {
        set_err(err, err_cap, "invalid image");
        return 0;
    }

    FILE* f = fopen(path, "wb");
    if (!f) {
        set_err(err, err_cap, "failed to open file for write");
        return 0;
    }

    unsigned char header[CHOOT_HEADER_SIZE];
    memset(header, 0, sizeof(header));
    memcpy(header, CHOOT_MAGIC, 8);
    header[8] = (unsigned char)(CHOOT_VERSION & 0xFF);
    header[9] = (unsigned char)((CHOOT_VERSION >> 8) & 0xFF);
    header[10] = 0;
    header[11] = 0;
    uint32_t atom_count = image->atom_count;
    header[12] = (unsigned char)(atom_count & 0xFF);
    header[13] = (unsigned char)((atom_count >> 8) & 0xFF);
    header[14] = (unsigned char)((atom_count >> 16) & 0xFF);
    header[15] = (unsigned char)((atom_count >> 24) & 0xFF);
    header[16] = (unsigned char)(CHOOT_HEADER_SIZE & 0xFF);
    header[17] = (unsigned char)((CHOOT_HEADER_SIZE >> 8) & 0xFF);
    header[18] = (unsigned char)((CHOOT_HEADER_SIZE >> 16) & 0xFF);
    header[19] = (unsigned char)((CHOOT_HEADER_SIZE >> 24) & 0xFF);

    if (fwrite(header, 1, CHOOT_HEADER_SIZE, f) != CHOOT_HEADER_SIZE) {
        fclose(f);
        set_err(err, err_cap, "failed to write header");
        return 0;
    }

    for (uint32_t i = 0; i < atom_count; ++i) {
        const ChootAtom* a = &image->atoms[i];
        unsigned char rec[20];
        uint16_t x = f32_to_f16(a->x);
        uint16_t y = f32_to_f16(a->y);
        uint16_t sxx = f32_to_f16(a->sxx);
        uint16_t sxy = f32_to_f16(a->sxy);
        uint16_t syy = f32_to_f16(a->syy);
        uint16_t alpha = f32_to_f16(a->alpha);
        uint16_t Y = f32_to_f16(a->Y);
        uint16_t Co = f32_to_f16(a->Co);
        uint16_t Cg = f32_to_f16(a->Cg);
        uint16_t flags = a->flags;

        rec[0] = (unsigned char)(x & 0xFF);
        rec[1] = (unsigned char)((x >> 8) & 0xFF);
        rec[2] = (unsigned char)(y & 0xFF);
        rec[3] = (unsigned char)((y >> 8) & 0xFF);
        rec[4] = (unsigned char)(sxx & 0xFF);
        rec[5] = (unsigned char)((sxx >> 8) & 0xFF);
        rec[6] = (unsigned char)(sxy & 0xFF);
        rec[7] = (unsigned char)((sxy >> 8) & 0xFF);
        rec[8] = (unsigned char)(syy & 0xFF);
        rec[9] = (unsigned char)((syy >> 8) & 0xFF);
        rec[10] = (unsigned char)(alpha & 0xFF);
        rec[11] = (unsigned char)((alpha >> 8) & 0xFF);
        rec[12] = (unsigned char)(Y & 0xFF);
        rec[13] = (unsigned char)((Y >> 8) & 0xFF);
        rec[14] = (unsigned char)(Co & 0xFF);
        rec[15] = (unsigned char)((Co >> 8) & 0xFF);
        rec[16] = (unsigned char)(Cg & 0xFF);
        rec[17] = (unsigned char)((Cg >> 8) & 0xFF);
        rec[18] = (unsigned char)(flags & 0xFF);
        rec[19] = (unsigned char)((flags >> 8) & 0xFF);

        if (fwrite(rec, 1, sizeof(rec), f) != sizeof(rec)) {
            fclose(f);
            set_err(err, err_cap, "failed to write atom");
            return 0;
        }
    }

    fclose(f);
    return 1;
}

void choot_free(ChootImage* image) {
    if (!image) return;
    free(image->atoms);
    image->atoms = NULL;
    image->atom_count = 0;
}
