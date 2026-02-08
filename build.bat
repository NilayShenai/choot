@echo off
setlocal

REM ---- CHOOT tools build script ----
REM Run from workspace root

REM Auto-detect cl.exe
where cl >nul 2>nul
if not errorlevel 1 goto :have_cl

if exist "C:\Program Files\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvarsall.bat" (
    call "C:\Program Files\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvarsall.bat" x64 >nul 2>nul
    goto :have_cl
)
if exist "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvarsall.bat" (
    call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvarsall.bat" x64 >nul 2>nul
    goto :have_cl
)
if exist "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvarsall.bat" (
    call "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvarsall.bat" x64 >nul 2>nul
    goto :have_cl
)
if exist "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvarsall.bat" (
    call "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvarsall.bat" x64 >nul 2>nul
    goto :have_cl
)
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
exit /b 1

:have_cl

echo [1/3] Building png2choot...
cl /nologo /O2 /W4 /openmp /D _CRT_SECURE_NO_WARNINGS /I src tools\png2choot.c src\choot.c src\choot_render.c src\wic_png.c /Fe:png2choot.exe /link windowscodecs.lib ole32.lib
if errorlevel 1 goto :fail

echo [2/3] Building choot_raster...
cl /nologo /O2 /W4 /openmp /D _CRT_SECURE_NO_WARNINGS /I src tools\choot_raster.c src\choot.c src\choot_render.c /Fe:choot_raster.exe
if errorlevel 1 goto :fail

echo [3/3] Building tests...
cl /nologo /O2 /W4 /openmp /D _CRT_SECURE_NO_WARNINGS /I src tests\choot_tests.c src\choot.c src\choot_render.c /Fe:choot_tests.exe
if errorlevel 1 goto :fail

echo.
echo === All tools built ===
echo   png2choot.exe     - Convert PNG to CHOOT
echo   choot_raster.exe  - Rasterize CHOOT to PPM
echo   choot_tests.exe   - Run tests
goto :end

:fail
echo.
echo === BUILD FAILED ===
exit /b 1

:end
endlocal
