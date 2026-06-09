# DCT — Distribution of Capacitance Timescales

Electrochemical impedance spectroscopy (EIS) analysis toolkit for surface-immobilized redox systems. Inverts EIS data into a distribution of electron-transfer rate constants g(k) using a Maxwell admittance ladder, Tikhonov regularization, and non-negative least squares (NNLS).

Initial application: EAB aptamer biosensors (concentration-dependent folding from g(k) peak shifts).

---

## Prerequisites

| Requirement | Version tested | Notes |
|---|---|---|
| Wolfram Mathematica | tested on 13.1 |
| VS Code | Any recent | Windows only for Wolfbook |
| [Wolfbook VS Code extension](https://marketplace.visualstudio.com/items?itemName=wolfbook.wolfbook) | 2.7.14+ | Provides the `.wb` notebook UI and Wolfram kernel integration |

---

## Installing Wolfbook

1. Open VS Code → Extensions → search **Wolfbook** → Install.
2. Open `EAB_DCT_analysis.wb` in VS Code. Wolfbook opens it as an interactive notebook.

---

## Setting up Claude Code / MCP (optional, for AI-assisted working)

If you want Claude Code to read and edit cells live, Wolfbook exposes an MCP server that Claude Code can connect to. This requires manual configuration.

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

### Step 3 — Reload

Close and reopen VS Code, then restart the Claude Code panel. In Claude code, you can type `/mcp` and you should see something like this:

![alt text](image.png)


---

## Repository layout

```
DCT/
├── EAB_DCT_analysis.wb     ← main notebook — start here
├── DCT.wl                  ← top-level package loader
├── core/
│   ├── DataImport.wl       — reads .txt EIS files (freq, Re[Z], Im[Z])
│   ├── DCTKernel.wl        — Maxwell ladder + Tikhonov + NNLS solver
│   ├── DCTSpectrum.wl      — batch fitting driver
│   ├── KKCheck.wl          — Lin-KK validity check (Schönleber 2014)
│   └── NNLSFit.wl          — non-negative least squares
├── analysis/
│   ├── AdmittanceStats.wl  — per-frequency variance, CV, power-law fit
│   ├── PeakAnalysis.wl     — peak-finding, spectrum integration
│   └── Sensitivity.wl      — lambda / freq-range / weight-power sweeps
├── plots/                  — all plot functions (dark background, white-on-dark)
├── data/
│   ├── drift_03252026/     — long-term drift data (E1/, E2/, E3/ subdirs)
│   └── concentration/      — ligand titration data (E1_178uM.txt style names)
└── legacy/                 — original .m scripts kept for reference
```

---

## Entry point: `EAB_DCT_analysis.wb`

Open this file in VS Code with Wolfbook installed. Run cells top-to-bottom. The notebook has five sections:

### Section 1 — Setup
- **Cell 1**: Clears all package contexts and reloads from disk. Run this first whenever you edit a `.wl` file.
- **Cell 2 (params)**: The single source of truth for all fit parameters. Edit here, re-run, then re-run downstream cells. **Do not set these variables in later cells** — they will silently overwrite Cell 2's values and cause hard-to-diagnose mismatches.

Set these upfront (the only things you need to change for a new dataset):
```wolfram
dataDir          = "...\\drift_03252026";  (* path to drift data — E1/E2/E3 subdirs *)
concentrationDir = "...\\spike";           (* path to concentration data *)
repFileIdx       = 1;                      (* index of representative file for KK/sweep plots *)
fileDecimation   = 1;                     (* 1 = all drift files, N = every Nth — use a large number while tuning to more quickly iterate *)
binsPerDecade    = 15;                     (* tau grid density (bins per decade of frequency)*)
tauMinFactor     = 0.1;                    (* grid padding factor below 1/(2π fMax) *)
tauMaxFactor     = 10.;                    (* grid padding factor above 1/(2π fMin) *)
```

These are filled in as you work through Sections 2–3 (the notebook prints the suggested values):
```wolfram
fMinHz      = 0.5;      (* ← from KK check, Section 2 *)
fMaxHz      = 1000.;    (* ← from CV(f) / power-law inflection, Section 3 *)
weightPower = -0.75;    (* ← from power-law fit, Section 3 *)
lambdaND    = 1*^-3;    (* ← from L-curve corner, Section 3 *)
```

Set after inspecting the g(k) spectra:
```wolfram
kPeakRanges = {{30, 100}, {100, 500}};  (* unbound / bound k windows — concentration section *)
```

### Section 2 — Data Validation (Lin-KK)
Runs a Kramers–Kronig check to identify the reliable frequency window. The output `{kkFMin, kkFMax}` gives a lower bound for `fMinHz`. Copy the printed values into Cell 2.

### Section 3 — Parameter Tuning (drift data)
Uses the long-term drift dataset to pin down `weightPower`, `fMinHz`/`fMaxHz`, and `lambdaND` before touching the science data.

Workflow:
1. Run the **admittance stats** cell → look at CV(f) plot. Find the contiguous window where CV < 0.1 — that's your valid frequency range.
2. Run the **variance power-law** cell → `weightPower` is printed directly.
3. Run the **lambda sweep** + **L-curve** cells → pick the λ at the corner of the L-curve.
4. Update Cell 2 with those three values, re-run Cell 2, then re-run the batch fit.

> **Gotcha — params cell re-run order:** If you run Cell 2 *after* a cell that sets `lambdaND` to a different value, Cell 2 wins. Always re-run Cell 2 last when updating parameters.

### Section 4 — Concentration Analysis
Fits DCT spectra for all electrodes across the ligand concentration series. Results used downstream for fraction-folded and biophysical model fitting.

File naming convention required: `E1_178uM.txt`, `E2_50uM.txt`, etc. The electrode prefix (`E1`/`E2`/`E3`) and the numeric concentration before `uM` are parsed for sorting and labelling.

### Section 5 — Export
Set `exportFigure = True` to write all plots as PDFs to `data/figures/`. Default is `False`.

---

## Data file format

Plain-text, three columns, space or tab delimited:
```
frequency_Hz    Re_Z_ohm    Im_Z_ohm
```
No header line. Im[Z] should be **negative** in the inductive-free range (standard EIS convention). The importer reads these with `DCTDataImport`ImportEIS`.

---

## Biophysical model (concentration section)

Two Langmuir-type models are fit to fraction-folded vs [ligand]:

| Model | Equation | Parameters |
|---|---|---|
| A — no NF | `ff = f0 + (1-f0)/(1 + 10^(logKD - log10 c))` | KD only |
| B — with NF | `ff = f0 + (1-NF)(1-f0)/(1 + 10^(logKD - log10 c))` | KD, NF |

`f0` is the baseline fraction folded at zero ligand. `NF` (non-folding fraction) is the fraction of aptamers that never fold regardless of ligand concentration. Model A inflates KD when NF > 0; Model B is the physically correct form.

---

## Running tests

```wolfram
Get["tests/TestDataImport.wl"]
Get["tests/TestNNLS.wl"]
Get["tests/TestDCTKernel.wl"]
```

Or run them from the terminal:
```
wolframscript -file tests/TestDataImport.wl
```
