/*
 * CHOOT Viewer — baremetal Win32 + OpenGL + Dear ImGui
 *
 * Opens .choot files directly, renders the Gaussian field to a
 * texture at arbitrary resolution, and displays it with zoom/pan.
 *
 * Build:  see build.bat
 * Usage:  choot_viewer.exe [file.choot]
 */

#define WIN32_LEAN_AND_MEAN
#define COBJMACROS
#include <windows.h>
#include <shellapi.h>
#include <commdlg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "choot.h"
#include "choot_render.h"

/* Dear ImGui (C++ compiled, called from C via wrapper) */
#include "imgui.h"
#include "imgui_impl_win32.h"
#include "imgui_impl_opengl3.h"
#include "imgui_impl_opengl3_loader.h"

/* Prefer discrete GPU on hybrid systems (NVIDIA/AMD) */
#if defined(_MSC_VER)
__declspec(dllexport) unsigned long NvOptimusEnablement = 0x00000001;
__declspec(dllexport) int AmdPowerXpressRequestHighPerformance = 1;
#endif

/* GL constants not in ImGui's minimal loader */
#ifndef GL_NEAREST
#define GL_NEAREST 0x2600
#endif
#ifndef GL_LINEAR
#define GL_LINEAR  0x2601
#endif
#ifndef GL_CLAMP_TO_EDGE
#define GL_CLAMP_TO_EDGE 0x812F
#endif

/* ------------------------------------------------------------------ */
/* App state                                                           */
/* ------------------------------------------------------------------ */

static struct {
    ChootImage  image;
    int         loaded;
    char        filepath[MAX_PATH];

    uint8_t*    front_pixels;
    uint8_t*    back_pixels;
    int         render_w;
    int         render_h;
    GLuint      texture;
    int         tex_valid;

    float       zoom;
    float       pan_x;
    float       pan_y;
    int         dragging;
    POINT       drag_start;
    float       drag_pan_x;
    float       drag_pan_y;

    int         target_w;
    int         target_h;

    /* async render thread */
    HANDLE      render_thread;
    HANDLE      render_event;
    volatile LONG render_request;
    volatile LONG render_ready;
    volatile LONG render_exit;
    volatile LONG render_inflight;
    int         pending_w;
    int         pending_h;
    CRITICAL_SECTION render_cs;

    int         needs_render;

    double      render_time_ms;
} app;

/* ------------------------------------------------------------------ */
/* Rendering                                                           */
/* ------------------------------------------------------------------ */

static void request_render(int w, int h) {
    if (!app.loaded) return;

    if (w <= 0) w = 512;
    if (h <= 0) h = 512;
    if (w > 4096) w = 4096;
    if (h > 4096) h = 4096;

    /* Skip if already rendered at this resolution and no pending change */
    if (app.tex_valid && app.render_w == w && app.render_h == h && !app.render_inflight) {
        app.needs_render = 0;
        return;
    }

    app.pending_w = w;
    app.pending_h = h;
    InterlockedExchange(&app.render_request, 1);
    SetEvent(app.render_event);
    app.needs_render = 0;
}

static DWORD WINAPI render_thread_proc(LPVOID param) {
    (void)param;
    for (;;) {
        WaitForSingleObject(app.render_event, INFINITE);

        if (app.render_exit) break;

        for (;;) {
            LONG req = InterlockedExchange(&app.render_request, 0);
            if (req <= 0) break;

            InterlockedExchange(&app.render_inflight, 1);

            int w = app.pending_w;
            int h = app.pending_h;

            uint8_t* buf = (uint8_t*)realloc(app.back_pixels, (size_t)w * h * 4);
            if (!buf) break;
            app.back_pixels = buf;

            LARGE_INTEGER freq, t0, t1;
            QueryPerformanceFrequency(&freq);
            QueryPerformanceCounter(&t0);
            choot_render_rgba8(&app.image, w, h, app.back_pixels);
            QueryPerformanceCounter(&t1);

            EnterCriticalSection(&app.render_cs);
            app.render_w = w;
            app.render_h = h;
            app.render_time_ms = (double)(t1.QuadPart - t0.QuadPart) / (double)freq.QuadPart * 1000.0;
            app.render_ready = 1;
            LeaveCriticalSection(&app.render_cs);

            InterlockedExchange(&app.render_inflight, 0);
        }
    }
    return 0;
}

static void load_choot_file(const char* path) {
    ChootImage img;
    char err[256];
    if (!choot_load(path, &img, err, sizeof(err))) {
        MessageBoxA(NULL, err, "CHOOT Load Error", MB_OK | MB_ICONERROR);
        return;
    }
    if (app.loaded) choot_free(&app.image);
    app.image = img;
    app.loaded = 1;
    strncpy(app.filepath, path, MAX_PATH - 1);
    app.filepath[MAX_PATH - 1] = '\0';
    app.zoom = 1.0f;
    app.pan_x = 0.0f;
    app.pan_y = 0.0f;
    app.needs_render = 1;
}

static void open_file_dialog(HWND hwnd) {
    char path[MAX_PATH] = {0};
    OPENFILENAMEA ofn;
    memset(&ofn, 0, sizeof(ofn));
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = hwnd;
    ofn.lpstrFilter = "CHOOT Files (*.choot)\0*.choot\0All Files (*.*)\0*.*\0";
    ofn.lpstrFile = path;
    ofn.nMaxFile = MAX_PATH;
    ofn.Flags = OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST;
    if (GetOpenFileNameA(&ofn)) {
        load_choot_file(path);
    }
}

/* ------------------------------------------------------------------ */
/* ImGui forward decl for WndProc                                      */
/* ------------------------------------------------------------------ */

extern IMGUI_IMPL_API LRESULT ImGui_ImplWin32_WndProcHandler(HWND, UINT, WPARAM, LPARAM);

static HWND g_hwnd;
static HDC  g_hdc;

static LRESULT CALLBACK wnd_proc(HWND hwnd, UINT msg, WPARAM wp, LPARAM lp) {
    if (ImGui_ImplWin32_WndProcHandler(hwnd, msg, wp, lp))
        return 1;

    switch (msg) {
    case WM_SIZE:
        return 0;
    case WM_DESTROY:
        PostQuitMessage(0);
        return 0;
    case WM_DROPFILES: {
        HDROP drop = (HDROP)wp;
        char path[MAX_PATH];
        if (DragQueryFileA(drop, 0, path, MAX_PATH)) {
            load_choot_file(path);
        }
        DragFinish(drop);
        return 0;
    }
    }
    return DefWindowProcA(hwnd, msg, wp, lp);
}

/* ------------------------------------------------------------------ */
/* WGL context creation                                                */
/* ------------------------------------------------------------------ */

static HGLRC create_gl_context(HDC hdc) {
    PIXELFORMATDESCRIPTOR pfd;
    memset(&pfd, 0, sizeof(pfd));
    pfd.nSize = sizeof(pfd);
    pfd.nVersion = 1;
    pfd.dwFlags = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER;
    pfd.iPixelType = PFD_TYPE_RGBA;
    pfd.cColorBits = 32;
    pfd.cDepthBits = 24;
    pfd.cStencilBits = 8;
    int fmt = ChoosePixelFormat(hdc, &pfd);
    SetPixelFormat(hdc, fmt, &pfd);
    HGLRC ctx = wglCreateContext(hdc);
    wglMakeCurrent(hdc, ctx);
    return ctx;
}

/* ------------------------------------------------------------------ */
/* Export PPM                                                          */
/* ------------------------------------------------------------------ */

static void export_ppm(void) {
    if (!app.loaded || !app.front_pixels) return;

    char path[MAX_PATH] = "output.ppm";
    OPENFILENAMEA ofn;
    memset(&ofn, 0, sizeof(ofn));
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = g_hwnd;
    ofn.lpstrFilter = "PPM Files (*.ppm)\0*.ppm\0All Files (*.*)\0*.*\0";
    ofn.lpstrFile = path;
    ofn.nMaxFile = MAX_PATH;
    ofn.Flags = OFN_OVERWRITEPROMPT;
    if (!GetSaveFileNameA(&ofn)) return;

    FILE* f = fopen(path, "wb");
    if (!f) return;
    fprintf(f, "P6\n%d %d\n255\n", app.render_w, app.render_h);
    for (int y = 0; y < app.render_h; ++y) {
        const uint8_t* row = app.front_pixels + (size_t)y * app.render_w * 4;
        for (int x = 0; x < app.render_w; ++x)
            fwrite(row + x * 4, 1, 3, f);
    }
    fclose(f);
}

/* ------------------------------------------------------------------ */
/* Main                                                                */
/* ------------------------------------------------------------------ */

int WINAPI WinMain(HINSTANCE inst, HINSTANCE prev, LPSTR cmdline, int show) {
    (void)prev;

    /* defaults */
    app.target_w = 512;
    app.target_h = 512;
    app.zoom = 1.0f;

    /* register window class */
    WNDCLASSEXA wc;
    memset(&wc, 0, sizeof(wc));
    wc.cbSize = sizeof(wc);
    wc.style = CS_OWNDC;
    wc.lpfnWndProc = wnd_proc;
    wc.hInstance = inst;
    wc.hCursor = LoadCursor(NULL, IDC_ARROW);
    wc.lpszClassName = "CHOOTViewer";
    RegisterClassExA(&wc);

    g_hwnd = CreateWindowExA(
        WS_EX_ACCEPTFILES,
        "CHOOTViewer", "CHOOT Viewer",
        WS_OVERLAPPEDWINDOW | WS_VISIBLE,
        CW_USEDEFAULT, CW_USEDEFAULT, 1100, 750,
        NULL, NULL, inst, NULL);

    g_hdc = GetDC(g_hwnd);
    HGLRC glctx = create_gl_context(g_hdc);

    /* enable vsync if available */
    {
        typedef BOOL (WINAPI *PFNWGLSWAPINTERVALEXTPROC)(int);
        PFNWGLSWAPINTERVALEXTPROC wglSwapIntervalEXT =
            (PFNWGLSWAPINTERVALEXTPROC)wglGetProcAddress("wglSwapIntervalEXT");
        if (wglSwapIntervalEXT) wglSwapIntervalEXT(1);
    }

    /* ImGui setup */
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    ImGui::StyleColorsDark();

    /* make it look good */
    ImGuiStyle& style = ImGui::GetStyle();
    style.WindowRounding = 6.0f;
    style.FrameRounding = 4.0f;
    style.GrabRounding = 4.0f;
    style.WindowPadding = ImVec2(10, 10);
    style.FramePadding = ImVec2(6, 4);
    style.ItemSpacing = ImVec2(8, 6);

    /* dark blue accent */
    ImVec4* colors = style.Colors;
    colors[ImGuiCol_WindowBg]         = ImVec4(0.08f, 0.08f, 0.12f, 1.0f);
    colors[ImGuiCol_TitleBg]          = ImVec4(0.06f, 0.06f, 0.10f, 1.0f);
    colors[ImGuiCol_TitleBgActive]    = ImVec4(0.12f, 0.14f, 0.24f, 1.0f);
    colors[ImGuiCol_FrameBg]          = ImVec4(0.12f, 0.12f, 0.18f, 1.0f);
    colors[ImGuiCol_FrameBgHovered]   = ImVec4(0.18f, 0.18f, 0.28f, 1.0f);
    colors[ImGuiCol_Button]           = ImVec4(0.20f, 0.22f, 0.36f, 1.0f);
    colors[ImGuiCol_ButtonHovered]    = ImVec4(0.28f, 0.30f, 0.48f, 1.0f);
    colors[ImGuiCol_ButtonActive]     = ImVec4(0.36f, 0.38f, 0.56f, 1.0f);
    colors[ImGuiCol_SliderGrab]       = ImVec4(0.40f, 0.45f, 0.70f, 1.0f);
    colors[ImGuiCol_SliderGrabActive] = ImVec4(0.50f, 0.55f, 0.80f, 1.0f);
    colors[ImGuiCol_Header]           = ImVec4(0.20f, 0.22f, 0.36f, 1.0f);
    colors[ImGuiCol_HeaderHovered]    = ImVec4(0.28f, 0.30f, 0.48f, 1.0f);
    colors[ImGuiCol_CheckMark]        = ImVec4(0.50f, 0.60f, 1.00f, 1.0f);
    colors[ImGuiCol_MenuBarBg]        = ImVec4(0.10f, 0.10f, 0.16f, 1.0f);

    ImGui_ImplWin32_InitForOpenGL(g_hwnd);
    ImGui_ImplOpenGL3_Init("#version 110");

    InitializeCriticalSection(&app.render_cs);
    app.render_exit = 0;
    app.render_inflight = 0;
    app.render_event = CreateEventA(NULL, FALSE, FALSE, NULL);
    app.render_thread = CreateThread(NULL, 0, render_thread_proc, NULL, 0, NULL);

    /* load file from command line if provided */
    if (cmdline && cmdline[0]) {
        char path[MAX_PATH];
        /* strip quotes */
        const char* src = cmdline;
        if (*src == '"') {
            src++;
            int i = 0;
            while (*src && *src != '"' && i < MAX_PATH - 1) path[i++] = *src++;
            path[i] = '\0';
        } else {
            strncpy(path, src, MAX_PATH - 1);
            path[MAX_PATH - 1] = '\0';
        }
        load_choot_file(path);
    }

    /* main loop */
    MSG msg;
    int running = 1;
    while (running) {
        while (PeekMessageA(&msg, NULL, 0, 0, PM_REMOVE)) {
            TranslateMessage(&msg);
            DispatchMessageA(&msg);
            if (msg.message == WM_QUIT) running = 0;
        }
        if (!running) break;

        /* render choot if needed */
        if (app.needs_render) request_render(app.target_w, app.target_h);

        /* if render finished, upload to texture */
        if (app.render_ready) {
            EnterCriticalSection(&app.render_cs);
            app.render_ready = 0;
            if (app.back_pixels) {
                uint8_t* tmp = app.front_pixels;
                app.front_pixels = app.back_pixels;
                app.back_pixels = tmp;
            }
            LeaveCriticalSection(&app.render_cs);

            if (app.front_pixels) {
                if (!app.tex_valid) {
                    glGenTextures(1, &app.texture);
                    app.tex_valid = 1;
                }
                glBindTexture(GL_TEXTURE_2D, app.texture);
                glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
                glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
                glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
                glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

                static int last_w = 0, last_h = 0;
                glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, app.render_w, app.render_h,
                             0, GL_RGBA, GL_UNSIGNED_BYTE, app.front_pixels);
                last_w = app.render_w;
                last_h = app.render_h;
                glBindTexture(GL_TEXTURE_2D, 0);
            }
        }

        /* begin ImGui frame */
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplWin32_NewFrame();
        ImGui::NewFrame();

        RECT cr;
        GetClientRect(g_hwnd, &cr);
        int win_w = cr.right - cr.left;
        int win_h = cr.bottom - cr.top;

        /* ---- Menu bar ---- */
        if (ImGui::BeginMainMenuBar()) {
            if (ImGui::BeginMenu("File")) {
                if (ImGui::MenuItem("Open...", "Ctrl+O")) open_file_dialog(g_hwnd);
                if (ImGui::MenuItem("Export PPM...", NULL, false, app.loaded)) export_ppm();
                ImGui::Separator();
                if (ImGui::MenuItem("Quit", "Alt+F4")) running = 0;
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("View")) {
                if (ImGui::MenuItem("Zoom to Fit", "F")) {
                    app.zoom = 1.0f;
                    app.pan_x = 0.0f;
                    app.pan_y = 0.0f;
                }
                if (ImGui::MenuItem("Zoom 100%", "1")) {
                    app.zoom = 1.0f;
                    app.pan_x = 0.0f;
                    app.pan_y = 0.0f;
                }
                if (ImGui::MenuItem("Zoom 200%", "2")) app.zoom = 2.0f;
                if (ImGui::MenuItem("Zoom 400%", "4")) app.zoom = 4.0f;
                ImGui::EndMenu();
            }
            ImGui::EndMainMenuBar();
        }

        /* ---- Controls panel ---- */
        ImGui::SetNextWindowPos(ImVec2(8, 30), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(280, 0), ImGuiCond_FirstUseEver);
        ImGui::Begin("Controls", NULL, ImGuiWindowFlags_AlwaysAutoResize);

        if (ImGui::Button("Open .choot", ImVec2(260, 32))) {
            open_file_dialog(g_hwnd);
        }

        ImGui::Separator();
        ImGui::Text("Render Resolution");
        int changed = 0;
        ImGui::SliderInt("Width",  &app.target_w, 32, 4096);
        int width_done  = ImGui::IsItemDeactivatedAfterEdit();
        ImGui::SliderInt("Height", &app.target_h, 32, 4096);
        int height_done = ImGui::IsItemDeactivatedAfterEdit();
        changed |= (width_done || height_done);

        /* preset buttons */
        if (ImGui::Button("256x256"))   { app.target_w = app.target_h = 256;  changed = 1; }
        ImGui::SameLine();
        if (ImGui::Button("512x512"))   { app.target_w = app.target_h = 512;  changed = 1; }
        ImGui::SameLine();
        if (ImGui::Button("1024x1024")) { app.target_w = app.target_h = 1024; changed = 1; }

        if (ImGui::Button("2048x2048")) { app.target_w = app.target_h = 2048; changed = 1; }
        ImGui::SameLine();
        if (ImGui::Button("4096x4096")) { app.target_w = app.target_h = 4096; changed = 1; }

        if (changed && app.loaded) app.needs_render = 1;

        if (ImGui::Button("Re-render", ImVec2(260, 28)) && app.loaded)
            app.needs_render = 1;

        ImGui::Separator();
        ImGui::Text("Zoom: %.1fx", app.zoom);
        ImGui::SliderFloat("##zoom", &app.zoom, 0.1f, 16.0f, "%.1f");

        if (ImGui::Button("Export PPM...", ImVec2(260, 28)) && app.loaded)
            export_ppm();

        ImGui::End();

        /* ---- Info panel ---- */
        if (app.loaded) {
            ImGui::SetNextWindowPos(ImVec2(8, (float)win_h - 110), ImGuiCond_Always);
            ImGui::SetNextWindowSize(ImVec2(280, 0), ImGuiCond_Always);
            ImGui::Begin("Info", NULL,
                ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
                ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoTitleBar);

            /* extract just filename */
            const char* fname = app.filepath;
            const char* p = strrchr(app.filepath, '\\');
            if (p) fname = p + 1;
            p = strrchr(fname, '/');
            if (p) fname = p + 1;

            ImGui::Text("File: %s", fname);
            ImGui::Text("Atoms: %u", app.image.atom_count);
            ImGui::Text("Rendered: %dx%d", app.render_w, app.render_h);
            ImGui::Text("Render time: %.1f ms", app.render_time_ms);
            ImGui::Text("File size: %u bytes",
                (unsigned)(CHOOT_HEADER_SIZE + app.image.atom_count * 20));
            ImGui::End();
        }

        /* ---- Image viewport ---- */
        if (app.loaded && app.tex_valid) {
            /* handle mouse wheel zoom on image area */
            if (!io.WantCaptureMouse) {
                if (io.MouseWheel != 0.0f) {
                    float factor = (io.MouseWheel > 0) ? 1.15f : (1.0f / 1.15f);
                    app.zoom *= factor;
                    if (app.zoom < 0.1f) app.zoom = 0.1f;
                    if (app.zoom > 16.0f) app.zoom = 16.0f;
                }
                /* pan with middle or right mouse button */
                if (ImGui::IsMouseDragging(ImGuiMouseButton_Right, 1.0f) ||
                    ImGui::IsMouseDragging(ImGuiMouseButton_Middle, 1.0f)) {
                    ImVec2 delta = io.MouseDelta;
                    app.pan_x += delta.x;
                    app.pan_y += delta.y;
                }
            }

            float img_w = app.render_w * app.zoom;
            float img_h = app.render_h * app.zoom;
            float panel_left = 296.0f; /* controls panel width + margin */
            float avail_w = (float)win_w - panel_left;
            float avail_h = (float)win_h - 30.0f; /* menu bar */

            float cx = panel_left + avail_w * 0.5f + app.pan_x;
            float cy = 30.0f + avail_h * 0.5f + app.pan_y;

            ImGui::GetBackgroundDrawList()->AddRectFilled(
                ImVec2(panel_left, 30.0f),
                ImVec2((float)win_w, (float)win_h),
                IM_COL32(20, 20, 28, 255));

            /* checkerboard bg */
            float ix0 = cx - img_w * 0.5f;
            float iy0 = cy - img_h * 0.5f;
            ImGui::GetBackgroundDrawList()->AddRectFilled(
                ImVec2(ix0, iy0), ImVec2(ix0 + img_w, iy0 + img_h),
                IM_COL32(40, 40, 48, 255));

            ImGui::GetBackgroundDrawList()->AddImage(
                (ImTextureID)(intptr_t)app.texture,
                ImVec2(ix0, iy0),
                ImVec2(ix0 + img_w, iy0 + img_h));

            /* border */
            ImGui::GetBackgroundDrawList()->AddRect(
                ImVec2(ix0 - 1, iy0 - 1),
                ImVec2(ix0 + img_w + 1, iy0 + img_h + 1),
                IM_COL32(80, 80, 120, 200));
        } else {
            /* no file loaded — show hint */
            float panel_left = 296.0f;
            float avail_w = (float)win_w - panel_left;
            float avail_h = (float)win_h - 30.0f;
            float cx = panel_left + avail_w * 0.5f;
            float cy = 30.0f + avail_h * 0.5f;

            ImGui::GetBackgroundDrawList()->AddRectFilled(
                ImVec2(panel_left, 30.0f),
                ImVec2((float)win_w, (float)win_h),
                IM_COL32(20, 20, 28, 255));

            const char* hint = "Drop a .choot file here or click Open";
            ImVec2 ts = ImGui::CalcTextSize(hint);
            ImGui::GetForegroundDrawList()->AddText(
                ImVec2(cx - ts.x * 0.5f, cy - ts.y * 0.5f),
                IM_COL32(100, 100, 140, 200), hint);
        }

        /* ---- keyboard shortcuts ---- */
        if (!io.WantCaptureKeyboard) {
            if (io.KeyCtrl && ImGui::IsKeyPressed(ImGuiKey_O))
                open_file_dialog(g_hwnd);
            if (ImGui::IsKeyPressed(ImGuiKey_F)) {
                app.zoom = 1.0f; app.pan_x = 0.0f; app.pan_y = 0.0f;
            }
            if (ImGui::IsKeyPressed(ImGuiKey_1)) app.zoom = 1.0f;
            if (ImGui::IsKeyPressed(ImGuiKey_2)) app.zoom = 2.0f;
            if (ImGui::IsKeyPressed(ImGuiKey_4)) app.zoom = 4.0f;
        }

        /* ---- GL render ---- */
        glViewport(0, 0, win_w, win_h);
        glClearColor(0.08f, 0.08f, 0.12f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        SwapBuffers(g_hdc);

        /* keep UI responsive; vsync handles pacing */
        Sleep(1);
    }

    /* cleanup */
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplWin32_Shutdown();
    ImGui::DestroyContext();

    wglMakeCurrent(NULL, NULL);
    wglDeleteContext(glctx);
    ReleaseDC(g_hwnd, g_hdc);

    if (app.loaded) choot_free(&app.image);
    app.render_exit = 1;
    if (app.render_event) SetEvent(app.render_event);
    if (app.render_thread) {
        WaitForSingleObject(app.render_thread, 2000);
        CloseHandle(app.render_thread);
    }
    if (app.render_event) CloseHandle(app.render_event);
    DeleteCriticalSection(&app.render_cs);
    free(app.front_pixels);
    free(app.back_pixels);

    return 0;
}
