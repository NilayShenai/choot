#include "wic_png.h"

#include <windows.h>
#include <wincodec.h>
#include <stdlib.h>
#include <string.h>

static void set_err(char* err, size_t cap, const char* msg) {
    if (!err || cap == 0) return;
    strncpy(err, msg, cap - 1);
    err[cap - 1] = '\0';
}

int wic_load_png(const char* path, WicImage* out_image, char* err, size_t err_cap) {
    if (!out_image) {
        set_err(err, err_cap, "null out_image");
        return 0;
    }
    memset(out_image, 0, sizeof(*out_image));

    HRESULT hr = CoInitializeEx(NULL, COINIT_MULTITHREADED);
    if (FAILED(hr)) {
        set_err(err, err_cap, "CoInitializeEx failed");
        return 0;
    }

    IWICImagingFactory* factory = NULL;
    hr = CoCreateInstance(&CLSID_WICImagingFactory, NULL, CLSCTX_INPROC_SERVER, &IID_IWICImagingFactory, (LPVOID*)&factory);
    if (FAILED(hr) || !factory) {
        CoUninitialize();
        set_err(err, err_cap, "WIC factory creation failed");
        return 0;
    }

    IWICBitmapDecoder* decoder = NULL;
    wchar_t wpath[MAX_PATH];
    int wlen = MultiByteToWideChar(CP_UTF8, 0, path, -1, wpath, MAX_PATH);
    if (wlen == 0) {
        factory->lpVtbl->Release(factory);
        CoUninitialize();
        set_err(err, err_cap, "path conversion failed");
        return 0;
    }

    hr = factory->lpVtbl->CreateDecoderFromFilename(factory, wpath, NULL, GENERIC_READ, WICDecodeMetadataCacheOnDemand, &decoder);
    if (FAILED(hr) || !decoder) {
        factory->lpVtbl->Release(factory);
        CoUninitialize();
        set_err(err, err_cap, "decoder creation failed");
        return 0;
    }

    IWICBitmapFrameDecode* frame = NULL;
    hr = decoder->lpVtbl->GetFrame(decoder, 0, &frame);
    if (FAILED(hr) || !frame) {
        decoder->lpVtbl->Release(decoder);
        factory->lpVtbl->Release(factory);
        CoUninitialize();
        set_err(err, err_cap, "frame decode failed");
        return 0;
    }

    UINT width = 0, height = 0;
    hr = frame->lpVtbl->GetSize(frame, &width, &height);
    if (FAILED(hr) || width == 0 || height == 0) {
        frame->lpVtbl->Release(frame);
        decoder->lpVtbl->Release(decoder);
        factory->lpVtbl->Release(factory);
        CoUninitialize();
        set_err(err, err_cap, "GetSize failed");
        return 0;
    }

    IWICFormatConverter* converter = NULL;
    hr = factory->lpVtbl->CreateFormatConverter(factory, &converter);
    if (FAILED(hr) || !converter) {
        frame->lpVtbl->Release(frame);
        decoder->lpVtbl->Release(decoder);
        factory->lpVtbl->Release(factory);
        CoUninitialize();
        set_err(err, err_cap, "CreateFormatConverter failed");
        return 0;
    }

    hr = converter->lpVtbl->Initialize(converter, (IWICBitmapSource*)frame, &GUID_WICPixelFormat32bppRGBA,
                                       WICBitmapDitherTypeNone, NULL, 0.0, WICBitmapPaletteTypeCustom);
    if (FAILED(hr)) {
        converter->lpVtbl->Release(converter);
        frame->lpVtbl->Release(frame);
        decoder->lpVtbl->Release(decoder);
        factory->lpVtbl->Release(factory);
        CoUninitialize();
        set_err(err, err_cap, "converter init failed");
        return 0;
    }

    size_t stride = (size_t)width * 4;
    size_t size = stride * (size_t)height;
    uint8_t* pixels = (uint8_t*)malloc(size);
    if (!pixels) {
        converter->lpVtbl->Release(converter);
        frame->lpVtbl->Release(frame);
        decoder->lpVtbl->Release(decoder);
        factory->lpVtbl->Release(factory);
        CoUninitialize();
        set_err(err, err_cap, "out of memory");
        return 0;
    }

    hr = converter->lpVtbl->CopyPixels(converter, NULL, (UINT)stride, (UINT)size, pixels);
    if (FAILED(hr)) {
        free(pixels);
        converter->lpVtbl->Release(converter);
        frame->lpVtbl->Release(frame);
        decoder->lpVtbl->Release(decoder);
        factory->lpVtbl->Release(factory);
        CoUninitialize();
        set_err(err, err_cap, "CopyPixels failed");
        return 0;
    }

    converter->lpVtbl->Release(converter);
    frame->lpVtbl->Release(frame);
    decoder->lpVtbl->Release(decoder);
    factory->lpVtbl->Release(factory);
    CoUninitialize();

    out_image->width = (uint32_t)width;
    out_image->height = (uint32_t)height;
    out_image->rgba = pixels;
    return 1;
}

void wic_free(WicImage* image) {
    if (!image) return;
    free(image->rgba);
    image->rgba = NULL;
    image->width = 0;
    image->height = 0;
}
