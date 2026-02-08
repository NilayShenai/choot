@echo off
setlocal

REM ---- CHOOT Viewer build script (MSVC) ----
REM Run from: viewer\ directory
REM Auto-detects Visual Studio and sets up cl.exe

REM Check if cl is already available
where cl >nul 2>nul
if not errorlevel 1 goto :have_cl

REM Try VS 2022 BuildTools
if exist "C:\Program Files\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvarsall.bat" (
    call "C:\Program Files\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvarsall.bat" x64 >nul 2>nul
    goto :have_cl
)
REM Try VS 2022 Community
if exist "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvarsall.bat" (
    call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvarsall.bat" x64 >nul 2>nul
    goto :have_cl
)
REM Try VS 2022 Professional
if exist "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvarsall.bat" (
    call "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvarsall.bat" x64 >nul 2>nul
    goto :have_cl
)
REM Try VS 2022 Enterprise
if exist "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvarsall.bat" (
    call "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvarsall.bat" x64 >nul 2>nul
    goto :have_cl
)
REM Try VS 2019
if exist "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvarsall.bat" (
    call "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvarsall.bat" x64 >nul 2>nul
    goto :have_cl
)
if exist "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvarsall.bat" (
    call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvarsall.bat" x64 >nul 2>nul
    goto :have_cl
)
REM Try VS 18 (Preview/Latest)
if exist "C:\Program Files (x86)\Microsoft Visual Studio\18\BuildTools\VC\Auxiliary\Build\vcvarsall.bat" (
    call "C:\Program Files (x86)\Microsoft Visual Studio\18\BuildTools\VC\Auxiliary\Build\vcvarsall.bat" x64 >nul 2>nul
    goto :have_cl
)

REM Last resort: search for vcvarsall.bat anywhere
for /f "delims=" %%i in ('dir /s /b "C:\Program Files*\vcvarsall.bat" 2^>nul') do (
    call "%%i" x64 >nul 2>nul
    goto :have_cl
)

echo ERROR: Could not find Visual Studio or Build Tools.
echo Install "Desktop development with C++" from:
echo   https://visualstudio.microsoft.com/visual-cpp-build-tools/
exit /b 1

:have_cl

if not exist imgui\imgui.cpp (
    echo ERROR: Dear ImGui not found. Run setup.ps1 first.
    echo   powershell -ExecutionPolicy Bypass -File setup.ps1
    exit /b 1
)

echo [1/3] Compiling ImGui...
cl /nologo /c /O2 /EHsc /W3 ^
    /I imgui ^
    /D _CRT_SECURE_NO_WARNINGS ^
    imgui\imgui.cpp ^
    imgui\imgui_draw.cpp ^
    imgui\imgui_tables.cpp ^
    imgui\imgui_widgets.cpp ^
    imgui\imgui_impl_win32.cpp ^
    imgui\imgui_impl_opengl3.cpp
if errorlevel 1 goto :fail

echo [2/3] Compiling CHOOT library...
cl /nologo /c /O2 /W4 /openmp ^
    /D _CRT_SECURE_NO_WARNINGS ^
    /I ..\src ^
    ..\src\choot.c ^
    ..\src\choot_render.c
if errorlevel 1 goto :fail

echo [3/3] Compiling and linking viewer...
cl /nologo /c /O2 /EHsc /W3 /TP ^
    /D _CRT_SECURE_NO_WARNINGS ^
    /I imgui ^
    /I ..\src ^
    viewer.c
if errorlevel 1 goto :fail

link /nologo /OUT:choot_viewer.exe ^
    viewer.obj ^
    imgui.obj imgui_draw.obj imgui_tables.obj imgui_widgets.obj ^
    imgui_impl_win32.obj imgui_impl_opengl3.obj ^
    choot.obj choot_render.obj ^
    user32.lib gdi32.lib opengl32.lib shell32.lib comdlg32.lib dwmapi.lib
if errorlevel 1 goto :fail

echo.
echo === Build OK: choot_viewer.exe ===
echo Usage: choot_viewer.exe [file.choot]
echo   or drag-and-drop a .choot file onto the window
goto :end

:fail
echo.
echo === BUILD FAILED ===
exit /b 1

:end
endlocal
