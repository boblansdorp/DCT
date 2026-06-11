# DCT — Distribution of Capacitance Timescales

Electrochemical impedance spectroscopy (EIS) analysis toolkit for surface-immobilized redox systems. Inverts EIS data into a distribution of electron-transfer rate constants g(k) using a Maxwell admittance ladder, Tikhonov regularization, and non-negative least squares (NNLS).

Initial application: EAB aptamer biosensors (concentration-dependent folding from g(k) peak shifts).

---

## Prerequisites

| Requirement | Version tested | Notes |
|---|---|---|
| Wolfram Mathematica | tested on 13.1 | |
| VS Code | Any recent | Windows only for Wolfbook |
| [Wolfbook VS Code extension](https://marketplace.visualstudio.com/items?itemName=wolfbook.wolfbook) | 2.7.14+ | Provides the `.wb` notebook UI and Wolfram kernel integration |
| [Node.js](https://nodejs.org/) | 18+ | Required for the Wolfbook MCP bridge (not installed by default on Windows) |

> **Note:** `.wb` notebooks require VS Code + Wolfbook. They cannot be opened in Mathematica directly. The `.wl` package files are standard Wolfram Language and work anywhere.

---

## Installing Wolfbook

1. Open VS Code → Extensions → search **Wolfbook** → Install.
2. Open `EAB_DCT_analysis.wb` in VS Code. Wolfbook opens it as an interactive notebook.

---

## Setting up Claude Code / MCP (optional, for AI-assisted working)

If you want Claude Code to read and edit cells live, Wolfbook exposes an MCP server that Claude Code can connect to.

### Step 1 — Find the bridge script path

After installing the extension, the bridge is at:
```
%USERPROFILE%\.vscode\extensions\wolfbook.wolfbook-<VERSION>-win32-x64\out\extension\claude-mcp\stdio-bridge.js
```
Replace `<VERSION>` with your installed version (check the extensions folder).

### Step 2 — Edit `~/.claude.json`

Add Wolfbook to **both** the global `mcpServers` block and the project-scoped block. It needs to be in both — the project scope alone is not picked up by the VSCode extension panel.

```json
{
  "mcpServers": {
    "wolfbook": {
      "type": "stdio",
      "command": "node",
      "args": ["C:\\Users\\<YOU>\\.vscode\\extensions\\wolfbook.wolfbook-2.7.14-win32-x64\\out\\extension\\claude-mcp\\stdio-bridge.js"],
      "env": {}
    }
  },
  "projects": {
    "C:\\Users\\<YOU>\\Documents\\DCT": {
      "mcpServers": {
        "wolfbook": {
          "type": "stdio",
          "command": "node",
          "args": ["C:\\Users\\<YOU>\\.vscode\\extensions\\wolfbook.wolfbook-2.7.14-win32-x64\\out\\extension\\claude-mcp\\stdio-bridge.js"],
          "env": {}
        }
      }
    }
  }
}
```

### Step 3 — Activate the MCP server in Wolfbook

With a `.wb` file open, press `Ctrl+Shift+P` and run **Wolfbook Claude MCP**. This starts the bridge process that Claude connects to.

### Step 4 — Reload Claude

Close and reopen VS Code, or restart the Claude Code panel. It may take **2–3 restart cycles** before the MCP server is recognised — this is a known Wolfbook quirk. In Claude Code, type `/mcp` to confirm the `wolfbook` server appears as connected.

![alt text](image.png)

### Known limitations

**Working directory differs from Mathematica.** `NotebookDirectory[]` is not available in Wolfbook — it returns `$Failed`. Use explicit absolute paths (as in cell 0.2) or `$InputFileName` only if the notebook was opened via `Get`. For path-relative operations, set `dataDir` explicitly.

---

## Repository layout

```
DCT/
├── EAB_DCT_analysis.wb     ← main analysis notebook — start here
├── runTests.wb             ← test suite notebook
├── DCT.wl                  ← top-level package loader
├── core/
│   ├── DataImport.wl       — reads .txt EIS files (freq, Re[Z], Im[Z])
│   ├── DCTKernel.wl        — Maxwell ladder + Tikhonov + NNLS solver
│   ├── DCTSpectrum.wl      — batch fitting driver
│   ├── KKCheck.wl          — Lin-KK validity check (Schönleber 2014)
│   └── NNLSFit.wl          — non-negative least squares
├── analysis/
│   ├── AdmittanceStats.wl  — per-frequency variance, CV, power-law fit
│   ├── Biophysics.wl       — 3-state/4-state aptamer folding model (symbolic + fitting)
│   ├── PeakAnalysis.wl     — peak-finding, spectrum integration
│   └── Sensitivity.wl      — lambda / freq-range / weight-power sweeps
├── plots/                  — all plot functions (theme-aware: dark + publication)
│   └── vendor/CustomTicks/ — vendored tick library, retained (MIT); no longer on the
│                              active path — ThemeLogTicks/ThemeLinTicks in PlotTheme.wl
│                              now generate outward ticks self-contained (see Acknowledgments)
├── tests/
│   ├── TestDataImport.wl
│   ├── TestNNLS.wl
│   ├── TestDCTKernel.wl
│   ├── TestDCTInversion.wl
│   └── TestPlotTicks.wl    — outward-tick regression tests (run via runTests.wb)
├── data/
│   ├── drift_03252026/     — long-term drift data (E1/, E2/, E3/ subdirs)
│   └── spike/              — ligand titration data (E1_178uM.txt style names)
└── legacy/                 — original .m scripts kept for reference
```

---

## Entry point: `EAB_DCT_analysis.wb`

Open in VS Code with Wolfbook installed. Run cells top-to-bottom within each section. The notebook has seven numbered sections:

### Section 0 — Setup
- **0.1 Load packages**: Clears all package contexts and reloads from disk. Run this first whenever you edit a `.wl` file.
- **0.2 Parameters**: Global fit parameters — the single source of truth for values used across sections 4 and 5. Edit here and re-run before running the batch fits.

Global parameters in 0.2:
```wolfram
dataDir          = "...\drift_03252026";  (* path to drift data *)
concentrationDir = "...\spike";           (* path to concentration data *)
lambdaND         = 1e-3;   (* Tikhonov regularisation — from L-curve, Section 3 *)
fMaxHz           = 1000.;  (* upper frequency cutoff — from admittance stats, Section 2 *)
weightPower      = -0.75;  (* admittance weight exponent — from power-law fit, Section 2 *)
binsPerDecade    = 15;     (* tau grid density *)
tauMinFactor     = 0.1;
tauMaxFactor     = 10.;
```

`fMinHz` is **not** set here — it is determined automatically in Section 1 from the Lin-KK population analysis and written by cell 1.3.

### Section 1 — KK Data Validation
Runs Lin-KK on **all** drift files to determine `fMinHz` from the population distribution. Section-specific parameters at the top of cell 1.1:
```wolfram
kkResidualThreshold = 0.02;  (* point-by-point |residual|/|Y| cutoff; try 0.05 if fMinHz is too high *)
kkFMinPercentile    = 50;    (* percentile of fMin distribution to use; 50 = median *)
```
Cell 1.2 shows a histogram of all files' KK lower bounds — inspect for outliers before accepting the result. Cell 1.3 writes `fMinHz`.

### Section 2 — Admittance Statistics
Computes per-frequency mean and variance of |Y(f)| across all drift files. Used to derive `fMaxHz` (from CV inflection) and `weightPower` (from variance power-law fit). Update 0.2 with the printed values.

### Section 3 — Parameter Tuning
L-curve and sweep plots on the representative drift file. Used to choose `lambdaND`. Workflow:
1. Run **3.1 Tau grid** and **3.2 Lambda sweep** to compute spectra across λ values.
2. Run **3a.1 Plot L-curve** to find the corner.
3. Set `lambdaND` in 0.2, re-run 0.2, then verify with **3a.3 Verify λ**.

### Section 4 — Drift Fit
Batch fits all drift files. `fileDecimation` is set at the top of cell 4.1 (default 1 = all files; increase to subsample during tuning). Results feed into the drift plot cells (4a.1–4a.7).

### Section 5 — Concentration Analysis
Batch fits all concentration files using the parameters locked in 0.2. Cells 5.1–5.4 produce g(k) overlays, area-normalised spectra, and fraction-folded vs [ligand]. Section 5a fits the biophysical model.

### Section 6 — Export
Set `exportFigure = True` to write all plots as PDFs to a `figures/` folder next to `dataDir`.

---

## Data file format

Plain-text, tab-delimited, with a header line:
```
Freq/Hz    Re(Z)/Ohm    -Im(Z)/Ohm
1000.0     120.3        -45.2
...
```
The importer (`DCTDataImport`ImportEIS`) normalises header strings and accepts common variants (`Re(Z)`, `Zre`, `Z'`, etc.). Im[Z] should be **negative** in the inductive-free range (standard EIS convention); the `-Im(Z)` column sign convention means the stored values are positive.

---

## Biophysical model (Section 5a)

Two mechanistic kinetic models are fit to fraction-folded vs [ligand] using symbolic derivation in `analysis/Biophysics.wl`:

| Model | States | Parameters |
|---|---|---|
| 3-state | U ↔ F ↔ B | K_S (structural switching), K_D (ligand dissociation) |
| 4-state | U ↔ F ↔ B + NF (non-folding) | K_S, K_D, NF |

Baseline fraction folded at zero ligand emerges from the model as K_S / (1 + K_S). The 4-state NF parameter absorbs the saturation shortfall, so its K_D is typically much smaller than the 3-state K_D.

Both model equations are printed symbolically at runtime before fitting.

---

## Running tests

Open `runTests.wb` in VS Code and run all cells, or run individual test files from the terminal:

```
wolframscript -file tests/TestDataImport.wl
wolframscript -file tests/TestNNLS.wl
wolframscript -file tests/TestDCTKernel.wl
wolframscript -file tests/TestDCTInversion.wl
```

All tests should report 0 failures. Plot/tick regression tests live in
`tests/TestPlotTicks.wl` and are run from `runTests.wb`.

---

## Acknowledgments / third-party code

- **CustomTicks** — log-aware tick generation, vendored under
  `plots/vendor/CustomTicks/`. © 2021 Mark A. Caprio, MIT License (full text in
  `plots/vendor/CustomTicks/LICENSE.md`). Part of the LevelScheme package. Retained
  for reference; the live plots now use native `ListLogLinearPlot`/`ListLogLogPlot`
  with self-contained outward-tick generators (`plots/PlotTheme.wl`), so CustomTicks
  is no longer loaded.
- **NNLS** (`core/NNLSFit.wl`) — adapted from **Michael Woodhams's** Mathematica
  implementation of the **Lawson–Hanson** active-set algorithm (Lawson & Hanson,
  *Solving Least Squares Problems*, Prentice-Hall 1974 / SIAM 1995), posted to
  comp.soft-sys.math.mathematica on 2 Oct 2003. Woodhams **placed the code in the
  public domain** and asks that his authorship be acknowledged.
  [Original thread](https://groups.google.com/g/comp.soft-sys.math.mathematica/c/cHFiKQ8ssaI)
  · [Stack Exchange repost](https://mathematica.stackexchange.com/questions/269727).
