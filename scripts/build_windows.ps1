param(
    [string]$Compiler = "gcc"
)

$compilerCommand = Get-Command $Compiler -ErrorAction SilentlyContinue
if (-not $compilerCommand) {
    $msysCompiler = "C:\msys64\ucrt64\bin\gcc.exe"
    if (Test-Path -LiteralPath $msysCompiler) {
        $Compiler = $msysCompiler
    }
    else {
        throw "gcc was not found. Install MSYS2 UCRT64 gcc or pass -Compiler <path>."
    }
}
else {
    $Compiler = $compilerCommand.Source
}

$root = Split-Path -Parent $PSScriptRoot
$src = Join-Path $root "src"
$out = Join-Path $root "build"
New-Item -ItemType Directory -Force -Path $out | Out-Null

$arguments = @(
    "-O3",
    "-shared",
    "-Wl,--export-all-symbols",
    "-o",
    (Join-Path $out "deepc_boxqp_optimized.dll"),
    (Join-Path $src "deepc_boxqp_optimized_ctypes_wrapper.c"),
    "-lm"
)

& $Compiler @arguments

if ($LASTEXITCODE -ne 0) {
    throw "C solver build failed. Install gcc (for example MSYS2 UCRT64 gcc) and retry."
}

Write-Host "Built $out\deepc_boxqp_optimized.dll"
