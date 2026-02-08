$ErrorActionPreference = "Stop"

$IMGUI_VERSION = "v1.91.8"
$IMGUI_URL = "https://raw.githubusercontent.com/ocornut/imgui/$IMGUI_VERSION"

$dir = Split-Path -Parent $MyInvocation.MyCommand.Path
$imguiDir = Join-Path $dir "imgui"

if (!(Test-Path $imguiDir)) {
    New-Item -ItemType Directory -Path $imguiDir | Out-Null
}

$files = @(
    "imgui.h",
    "imgui.cpp",
    "imgui_draw.cpp",
    "imgui_tables.cpp",
    "imgui_widgets.cpp",
    "imgui_internal.h",
    "imconfig.h",
    "imstb_rectpack.h",
    "imstb_textedit.h",
    "imstb_truetype.h",
    "backends/imgui_impl_win32.h",
    "backends/imgui_impl_win32.cpp",
    "backends/imgui_impl_opengl3.h",
    "backends/imgui_impl_opengl3.cpp",
    "backends/imgui_impl_opengl3_loader.h"
)

foreach ($f in $files) {
    $url = "$IMGUI_URL/$f"
    $name = Split-Path -Leaf $f
    $dest = Join-Path $imguiDir $name
    if (Test-Path $dest) {
        Write-Host "  skip  $name (exists)"
        continue
    }
    Write-Host "  fetch $name"
    Invoke-WebRequest -Uri $url -OutFile $dest -UseBasicParsing
}

Write-Host ""
Write-Host "Dear ImGui $IMGUI_VERSION downloaded to viewer/imgui/"
Write-Host "Run build.bat to compile the viewer."
