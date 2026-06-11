(* ::Package:: *)

(*  DCT.wl  -  master loader
    Loads all DCT sub-packages from the modular directory structure.

    Primary entry point (unchanged from previous versions):
      DCT`DCTSpectrum[file]           -- fit a single EIS file
    Sub-package layout:
      core/DataImport.wl        EIS file parsing
      core/NNLSFit.wl           Non-negative least squares solver
      core/DCTKernel.wl         Maxwell ladder + Tikhonov kernel
      core/DCTSpectrum.wl       Outer Brent loop, DCTSpectrum[] entry point
      analysis/PeakAnalysis.wl  Peak picking, spectrum integration
      analysis/Sensitivity.wl   Parameter sweep utilities
      analysis/Biophysics.wl    3/4-state aptamer folding model (symbolic derivation + fitting)
      plots/PlotSpectrum.wl     g(k), g(tau), cumulative capacitance
      plots/PlotImpedance.wl    Nyquist, Bode magnitude/phase
      plots/PlotResiduals.wl    Fit residual analysis
      plots/PlotTimeSeries.wl   Rs, C0, kPeak vs time
      plots/PlotHeatmap.wl      2-D g(k) heatmap over time
    Tests (run independently):
      tests/TestDataImport.wl
      tests/TestNNLS.wl
      tests/TestDCTKernel.wl
      tests/TestPlots.wl
*)

Module[{here, load},

  here = DirectoryName[$InputFileName];
  load = Function[rel, Get[FileNameJoin[{here, rel}]]];

  (* Core: order matters — NNLSFit before Kernel, Kernel before Spectrum *)
  load["core/NNLSFit.wl"];
  load["core/DataImport.wl"];
  load["core/KKCheck.wl"];
  load["core/DCTKernel.wl"];
  load["core/DCTSpectrum.wl"];

  (* Analysis *)
  load["analysis/PeakAnalysis.wl"];
  load["analysis/Sensitivity.wl"];
  load["analysis/AdmittanceStats.wl"];
  load["analysis/Biophysics.wl"];

  (* Plots *)
  (* Mark Caprio, MIT; outward ticks. Quiet the shadow warnings: CustomTicks
     exports LogPlot/TickDirection (also System` symbols) — we always reference
     them as CustomTicks`LogPlot / CustomTicks`TickDirection. *)
  Quiet[load["plots/vendor/CustomTicks/CustomTicks.m"], {General::shdw}];
  load["plots/PlotTheme.wl"];
  load["plots/PlotSpectrum.wl"];
  load["plots/PlotImpedance.wl"];
  load["plots/PlotResiduals.wl"];
  load["plots/PlotTimeSeries.wl"];
  load["plots/PlotHeatmap.wl"];
  load["plots/PlotKK.wl"];
  load["plots/PlotAdmittanceStats.wl"];
  load["plots/PlotTuning.wl"];

  Print[Style["\:2714 DCT package loaded.", Darker[Green], 13]];
]
