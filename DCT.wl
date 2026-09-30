(* ::Package:: *)

(*  DCT.wl  -  main loader
    Loads all DCT sub-packages from the directory this file is in.

    Primary entry point:
      DCT`DCTSpectrum[file]           -- fit a single EIS file
    Sub-package layout:
      core/NNLSFit.wl                 Non-negative least squares solver (Lawson-Hanson)
      core/DataImport.wl              EIS file parsing
      core/KKCheck.wl                 Linearized Kramers-Kronig check in admittance space
      core/DCTKernel.wl               Maxwell ladder + Tikhonov kernel
      core/DCTSpectrum.wl             Outer Brent loop, DCTSpectrum[] entry point
      core/CPEFit.wl                  Equivalent-circuit (CPE) fit: Rs + Cdl||(Rct+CPE)
      analysis/PeakAnalysis.wl        Peak picking, spectrum integration
      analysis/Sensitivity.wl         Parameter sweep utilities
      analysis/AdmittanceStats.wl     Per-frequency admittance variance, power-law fit
      analysis/Biophysics.wl          3/4-state aptamer folding model (symbolic derivation + fitting)
      plots/PlotTheme.wl              Dark/publication themes, outward tick generators
      plots/PlotSpectrum.wl           g(k), g(tau), cumulative capacitance
      plots/PlotImpedance.wl          Nyquist, Bode magnitude/phase
      plots/PlotResiduals.wl          Fit residual analysis
      plots/PlotTimeSeries.wl         Rs, C0, kPeak vs time
      plots/PlotHeatmap.wl            2-D g(k) heatmap over time
      plots/PlotKK.wl                 Kramers-Kronig residuals and valid-range histogram
      plots/PlotAdmittanceStats.wl    Admittance variance vs mean
      plots/PlotTuning.wl             Lambda sweep and L-curve
      plots/PlotValidation.wl         Synthetic-recovery validation
    Tests (run from runTests.wb or with wolframscript -file):
      tests/TestNNLS.wl, TestDataImport.wl, TestDCTKernel.wl, TestDCTInversion.wl,
      tests/TestCPEFit.wl, TestBiophysics.wl, TestPlotTicks.wl
*)

Module[{here, load},

  here = DirectoryName[$InputFileName];
  load = Function[rel, Get[FileNameJoin[{here, rel}]]];

  (* Core: order matters: NNLSFit before Kernel, Kernel before Spectrum *)
  load["core/NNLSFit.wl"];
  load["core/DataImport.wl"];
  load["core/KKCheck.wl"];
  load["core/DCTKernel.wl"];
  load["core/DCTSpectrum.wl"];
  load["core/CPEFit.wl"];

  (* Analysis *)
  load["analysis/PeakAnalysis.wl"];
  load["analysis/Sensitivity.wl"];
  load["analysis/AdmittanceStats.wl"];
  load["analysis/Biophysics.wl"];

  (* Plots *)
  load["plots/PlotTheme.wl"];
  load["plots/PlotSpectrum.wl"];
  load["plots/PlotImpedance.wl"];
  load["plots/PlotResiduals.wl"];
  load["plots/PlotTimeSeries.wl"];
  load["plots/PlotHeatmap.wl"];
  load["plots/PlotKK.wl"];
  load["plots/PlotAdmittanceStats.wl"];
  load["plots/PlotTuning.wl"];
  load["plots/PlotValidation.wl"];

  Print[Style["\:2714 DCT package loaded.", Darker[Green], 13]];
]
