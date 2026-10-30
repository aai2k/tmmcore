<#
.SYNOPSIS
    Writes the Essential Macleod results that tests/macleod.mjs checks tmmcore
    against, one JSON file per case beside this script.

.DESCRIPTION
    Drives Essential Macleod's COM automation server (ProgID EMacleod.Session).

    Every material is a constant n and k defined below, written to Essential
    Macleod as a two-point table at 300 and 2500 nm under a name starting with
    "tmmcore", so both programs use the same optical constants at every
    wavelength. The values are binary fractions, exact in single and double
    precision. Each case is built in memory and calculated; no design file is
    saved. Re-running the script replaces the materials and the JSON files.

    For each case and angle of incidence, at every wavelength of the case's
    grid, in s and p: R, T, A, the reflectance seen from the substrate side,
    the phase of r and t, GD, GDD and TOD of r and t, and Psi and Delta of r and
    t. For the loci listed with a case, the reflectance through the stack at
    one wavelength, from the substrate interface to the outer surface: the
    monitoring curve of every layer in turn.

    Values are written as Essential Macleod reports them: indices as n - ik
    with k >= 0, R, T and A as fractions, phase, Psi and Delta in degrees,
    Delta in the Azzam-Bashara sign, GD in fs, GDD in fs^2, TOD in fs^3, depth
    in nm from the outer surface, the locus reflectance in percent. The layers
    are listed from the incident medium to the substrate, as tmmcore takes
    them. tests/macleod.mjs converts to tmmcore's conventions.

    Requires Essential Macleod with its automation server registered.

.PARAMETER OutDir
    Where the JSON files go. Defaults to this script's folder.
#>
param(
    [string]$OutDir = $PSScriptRoot
)

$ErrorActionPreference = 'Stop'
$ci = [System.Globalization.CultureInfo]::InvariantCulture

# Round-trip formatting keeps all 17 digits; JSON has no NaN, so a value
# Essential Macleod leaves undefined is written as null.
function Num([double]$v) {
    if ([double]::IsNaN($v) -or [double]::IsInfinity($v)) { return 'null' }
    return $v.ToString('R', $ci)
}
function NumArray($values) { return '[' + ((@($values) | ForEach-Object { Num $_ }) -join ',') + ']' }
function Str([string]$s) { return '"' + $s.Replace('\', '\\').Replace('"', '\"') + '"' }

# n and k of every material, k >= 0 for absorption.
$Materials = [ordered]@{
    'tmmcore Air'       = @(1.0, 0.0)
    'tmmcore L'         = @(1.4375, 0.0)
    'tmmcore H'         = @(2.25, 0.0)
    'tmmcore S'         = @(1.5, 0.0)
    'tmmcore K'         = @(2.0, 0.0625)
    'tmmcore M'         = @(0.125, 3.5)
    'tmmcore N0 weak'   = @(1.5, 0.0009765625)
    'tmmcore N0 strong' = @(1.5, 0.015625)
}

$Visible = @{ Min = 400.0; Max = 800.0; Step = 10.0 }
$Stack9 = @('tmmcore H', 62.25, 'tmmcore L', 97.5, 'tmmcore H', 62.25, 'tmmcore L', 97.5, 'tmmcore H', 62.25,
            'tmmcore L', 97.5, 'tmmcore H', 62.25, 'tmmcore L', 97.5, 'tmmcore H', 62.25)
$Stack3 = @('tmmcore H', 62.25, 'tmmcore L', 97.5, 'tmmcore H', 62.25)

# Layers as material, thickness in nm, pairs, from the incident medium.
$Cases = @(
    @{ Id = 'dielectric-stack'; Incident = 'tmmcore Air'; Substrate = 'tmmcore S'; Layers = $Stack9
       Angles = @(0.0, 30.0, 60.0); Grid = $Visible
       Loci = @(@{ Wavelength = 560.0; Angle = 0.0 }, @{ Wavelength = 450.0; Angle = 45.0 }) },
    @{ Id = 'metal-mirror'; Incident = 'tmmcore Air'; Substrate = 'tmmcore S'; Layers = @('tmmcore L', 82.5, 'tmmcore M', 120.0)
       Angles = @(0.0, 45.0, 70.0); Grid = $Visible
       Loci = @(@{ Wavelength = 560.0; Angle = 45.0 }) },
    @{ Id = 'absorbing-layers'; Incident = 'tmmcore Air'; Substrate = 'tmmcore S'
       Layers = @('tmmcore K', 80.0, 'tmmcore L', 100.0, 'tmmcore K', 80.0, 'tmmcore L', 100.0)
       Angles = @(0.0, 30.0, 60.0); Grid = $Visible
       Loci = @(@{ Wavelength = 500.0; Angle = 30.0 }) },
    @{ Id = 'absorbing-substrate'; Incident = 'tmmcore Air'; Substrate = 'tmmcore K'; Layers = @('tmmcore L', 97.5, 'tmmcore H', 62.25)
       Angles = @(0.0, 30.0, 60.0); Grid = $Visible
       Loci = @(@{ Wavelength = 600.0; Angle = 30.0 }) },
    @{ Id = 'absorbing-incident-weak'; Incident = 'tmmcore N0 weak'; Substrate = 'tmmcore S'; Layers = $Stack3
       Angles = @(0.0, 30.0, 60.0); Grid = $Visible; Loci = @() },
    @{ Id = 'absorbing-incident-strong'; Incident = 'tmmcore N0 strong'; Substrate = 'tmmcore S'; Layers = $Stack3
       Angles = @(0.0, 30.0, 60.0); Grid = $Visible; Loci = @() },
    @{ Id = 'frustrated-tir'; Incident = 'tmmcore S'; Substrate = 'tmmcore S'; Layers = @('tmmcore Air', 400.0)
       Angles = @(30.0, 45.0, 60.0); Grid = $Visible; Loci = @() },
    @{ Id = 'grazing'; Incident = 'tmmcore Air'; Substrate = 'tmmcore S'; Layers = @('tmmcore L', 97.5)
       Angles = @(80.0, 89.0, 89.9); Grid = $Visible; Loci = @() },
    @{ Id = 'opaque-metal'; Incident = 'tmmcore Air'; Substrate = 'tmmcore S'; Layers = @('tmmcore M', 1000.0)
       Angles = @(0.0, 45.0); Grid = $Visible
       Loci = @(@{ Wavelength = 560.0; Angle = 0.0 }) },
    @{ Id = 'thick-slab'; Incident = 'tmmcore Air'; Substrate = 'tmmcore S'; Layers = @('tmmcore L', 20000.0)
       Angles = @(0.0, 45.0); Grid = @{ Min = 500.0; Max = 600.0; Step = 2.5 }; Loci = @() }
)

# DesignXAxisType code of a wavelength scan.
$WavelengthAxis = 87

# DesignYAxisType codes and the key each is written under.
$YTypes = @(
    @(82, 'R'), @(84, 'T'), @(65, 'A'), @(66, 'Rback'),
    @(80, 'rPhaseDeg'), @(81, 'tPhaseDeg'),
    @(71, 'rGd'), @(72, 'rGdd'), @(73, 'rTod'), @(74, 'tGd'), @(75, 'tGdd'), @(76, 'tTod'),
    @(88, 'rPsiDeg'), @(87, 'rDeltaDeg'), @(86, 'tPsiDeg'), @(85, 'tDeltaDeg')
)

if (-not (Test-Path $OutDir)) { New-Item -ItemType Directory -Force -Path $OutDir | Out-Null }
# The program version, from the executable its COM server is registered as. It
# is a 32-bit server, so on 64-bit Windows its class sits under WOW6432Node.
$clsid = (Get-ItemProperty -LiteralPath 'Registry::HKEY_CLASSES_ROOT\EMacleod.Session\CLSID').'(default)'
$server = @("Registry::HKEY_CLASSES_ROOT\CLSID\$clsid\LocalServer32",
            "Registry::HKEY_CLASSES_ROOT\WOW6432Node\CLSID\$clsid\LocalServer32") |
    Where-Object { Test-Path -LiteralPath $_ } | Select-Object -First 1
$command = (Get-ItemProperty -LiteralPath $server).'(default)'
$exe = if ($command -match '^"([^"]+)"') { $Matches[1] } else { $command }
$version = (Get-Item -LiteralPath $exe).VersionInfo.FileVersion

Write-Host 'Opening COM session...'
$session = New-Object -ComObject EMacleod.Session
$library = $session.GetMaterials()

# ── Materials ────────────────────────────────────────────────────────────────

$existing = @($library.GetMaterialNames())
$used = @{}
foreach ($name in $Materials.Keys) {
    $nk = $Materials[$name]
    # Overwritten rather than deleted and recreated: a delete leaves an entry
    # behind that the next listing still shows.
    $m = if ($existing -contains $name) { $library.OpenMaterial($name, $false) } else { $library.NewMaterial() }
    $m.DeleteAll()
    $m.AddOpticalConstant(300.0, [double]$nk[0], [double]$nk[1])
    $m.AddOpticalConstant(2500.0, [double]$nk[0], [double]$nk[1])
    $m.SaveAs($name, $true)
    $back = $library.OpenMaterial($name, $true)
    $n = 0.0; $k = 0.0
    $back.GetOpticalConstant(550.0, [ref]$n, [ref]$k)
    if ($n -ne $nk[0] -or $k -ne $nk[1]) { throw "$name reads back as $n, $k" }
    $used[$name] = @($n, $k)
    Write-Host ("  {0,-18} n {1}  k {2}" -f $name, $n, $k)
}

# ── Cases ────────────────────────────────────────────────────────────────────

foreach ($case in $Cases) {
    $design = $session.NewDesign()
    $design.MediumName = $case.Incident
    $design.SubstrateName = $case.Substrate
    # A new design arrives with one quarterwave layer in it.
    if ($design.NumberOfLayers -gt 0) { $design.DeleteLayers(1, $design.NumberOfLayers) }
    $count = $case.Layers.Count / 2
    for ($i = 1; $i -le $count; $i++) { $design.InsertLayer($i) }
    for ($i = 1; $i -le $count; $i++) {
        $design.SetMaterialName($i, $case.Layers[2 * $i - 2])
        $design.SetPhysicalThickness($i, [double]$case.Layers[2 * $i - 1])
    }

    $layerJson = @()
    for ($i = 1; $i -le $count; $i++) {
        $name = ''; $d = 0.0
        $design.GetMaterialName($i, [ref]$name)
        $design.GetPhysicalThickness($i, [ref]$d)
        if ($d -ne [double]$case.Layers[2 * $i - 1]) { Write-Host "  $($case.Id) layer $i reads back as $d nm; that value is written" }
        $layerJson += '{"material":' + (Str $name) + ',"n":' + (NumArray $used[$name]) + ',"d":' + (Num $d) + '}'
    }

    $x = $design.XAxisParameters
    $y = $design.YAxisParameters
    if ($x.ConeWidth -ne 0) { throw "$($case.Id): cone width $($x.ConeWidth), not a collimated beam" }

    $scanJson = @()
    foreach ($angle in $case.Angles) {
        $x.Axis = $WavelengthAxis
        $x.Min = [double]$case.Grid.Min
        $x.Max = [double]$case.Grid.Max
        $x.TableInterval = [double]$case.Grid.Step
        $x.CalculationBandwidth = 0.0
        $x.CalculationIncidentAngle = [double]$angle
        $lambdas = $null
        $s = [ordered]@{}; $p = [ordered]@{}
        foreach ($t in $YTypes) {
            $y.Axis = $t[0]
            $ds = $design.CalculateTable($x, $y)
            $count2 = $ds.NumberOfItems
            if ($null -eq $lambdas) {
                $lambdas = New-Object 'double[]' $count2
                for ($i = 1; $i -le $count2; $i++) { $lambdas[$i - 1] = $ds.X($i) }
            }
            $vs = New-Object 'double[]' $count2; $vp = New-Object 'double[]' $count2
            for ($i = 1; $i -le $count2; $i++) { $vs[$i - 1] = $ds.Ys($i); $vp[$i - 1] = $ds.Yp($i) }
            $s[$t[1]] = $vs; $p[$t[1]] = $vp
        }
        $sJson = ($s.Keys | ForEach-Object { (Str $_) + ':' + (NumArray $s[$_]) }) -join ','
        $pJson = ($p.Keys | ForEach-Object { (Str $_) + ':' + (NumArray $p[$_]) }) -join ','
        $scanJson += '{"angle":' + (Num $angle) + ',"wavelengths":' + (NumArray $lambdas) + ',"s":{' + $sJson + '},"p":{' + $pJson + '}}'
    }

    $lociJson = @()
    $an = $design.Analysis
    foreach ($locus in $case.Loci) {
        foreach ($pol in @('s', 'p')) {
            $an.Wavelength = [double]$locus.Wavelength
            $an.IncidentAngle = [double]$locus.Angle
            $an.Polarization = if ($pol -eq 'p') { 1 } else { 2 }
            $a = $an.CalculateAdmittance()
            $depth = New-Object 'double[]' $a.NumberOfItems
            $rpct = New-Object 'double[]' $a.NumberOfItems
            for ($i = 1; $i -le $a.NumberOfItems; $i++) { $depth[$i - 1] = $a.X($i); $rpct[$i - 1] = $a.Reflectance($i) }
            $lociJson += '{"wavelength":' + (Num $locus.Wavelength) + ',"angle":' + (Num $locus.Angle) + ',"pol":' + (Str $pol) +
                ',"depth":' + (NumArray $depth) + ',"Rpct":' + (NumArray $rpct) + '}'
        }
    }

    $json = '{"case":' + (Str $case.Id) + ',"program":' + (Str "Essential Macleod $version") +
        ',"incident":{"material":' + (Str $case.Incident) + ',"n":' + (NumArray $used[$case.Incident]) + '}' +
        ',"substrate":{"material":' + (Str $case.Substrate) + ',"n":' + (NumArray $used[$case.Substrate]) + '}' +
        ",`n""layers"":[" + ($layerJson -join ',') + ']' +
        ",`n""scans"":[`n" + ($scanJson -join ",`n") + "`n]" +
        ",`n""loci"":[`n" + ($lociJson -join ",`n") + "`n]}`n"
    $path = Join-Path $OutDir ($case.Id + '.json')
    [System.IO.File]::WriteAllText($path, $json, (New-Object System.Text.UTF8Encoding $false))
    Write-Host ("  {0}: {1} layers, {2} angles, {3} loci -> {4}" -f $case.Id, $count, $case.Angles.Count, $lociJson.Count, $path)
}
