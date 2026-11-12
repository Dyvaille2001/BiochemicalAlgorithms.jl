# Run from any working directory; logs go under test/nma_best_final/logs.
$DataRoot = Split-Path -Parent $PSScriptRoot
$ProjectRoot = (Resolve-Path (Join-Path $DataRoot '..\..')).Path
$CodeDir = Join-Path $DataRoot 'code'
$LogDir = Join-Path $DataRoot 'logs'
New-Item -ItemType Directory -Force -Path $LogDir | Out-Null

function Step([string]$title, [scriptblock]$run, [string]$log) {
    Write-Host "`n=== $title ===" -ForegroundColor Cyan
    Push-Location $ProjectRoot
    try {
        & $run 2>&1 | Tee-Object -FilePath $log
        if ($LASTEXITCODE -ne 0) {
            Write-Host "    -> failed (exit code $LASTEXITCODE); see $log" -ForegroundColor Yellow
        } else {
            Write-Host "    -> OK, log: $log" -ForegroundColor Green
        }
    } finally {
        Pop-Location
    }
}

Step '1  NMA correctness checks' {
    julia --project=. (Join-Path $CodeDir 'test_nma_correctness.jl')
} (Join-Path $LogDir 'log_1_tests.txt')

Step '2  Amber force field vs MMTK references' {
    julia --project=. (Join-Path $CodeDir 'validate_forcefield.jl')
} (Join-Path $LogDir 'log_2_forcefield.txt')

Step '3  Lysozyme (1AKI) model validation' {
    julia --project=. (Join-Path $CodeDir 'validate_model.jl')
} (Join-Path $LogDir 'log_3_model.txt')

Write-Host "`nDone. See validation_forcefield, validation_model, and 1aki for results." -ForegroundColor Cyan
