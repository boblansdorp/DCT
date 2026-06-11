(* ::Package:: *)

(* DCTPlots`PlotTimeSeries
   Time-series plots of fit parameters extracted from a batch of spectra.
   x-axis is experiment time in hours (from spec["FinishTimeS"]).
   Each electrode subdirectory is shown in a distinct colour.

   Public API:
     PlotRsVsTime[goodSpecs]      solution resistance Rs vs time
     PlotC0VsTime[goodSpecs]      fast-limit capacitance C0 vs time
     PlotKPeakVsTime[goodSpecs]   dominant peak rate k vs time
     PlotTimeSeries[goodSpecs]    GraphicsGrid with all three panels

   Theme-aware: default DCTPlots`$DCTTheme ("Dark" | "Publication"); override
   one call with Theme -> "Publication". Electrode colours are theme-independent.
*)

BeginPackage["DCTPlots`", {"DCTPeakAnalysis`"}]

PlotRsVsTime::usage =
  "PlotRsVsTime[goodSpecs] plots Rs (\[CapitalOmega]) vs time (hours), \
one trace per electrode."

PlotC0VsTime::usage =
  "PlotC0VsTime[goodSpecs] plots C0 (nF) vs time (hours), \
one trace per electrode."

PlotKPeakVsTime::usage =
  "PlotKPeakVsTime[goodSpecs] plots the global dominant peak k (s^-1) vs time. \
PlotKPeakVsTime[goodSpecs, {{kMin1,kMax1},{kMin2,kMax2},...}] shows one panel \
per k-range, tracking the peak within each window."

PlotTimeSeries::usage =
  "PlotTimeSeries[goodSpecs] returns a GraphicsGrid with Rs, C0, and kPeak \
vs time panels, each electrode in a distinct colour."

Begin["`Private`"]

tsImageSize = 380;
tsAspect    = 0.75;

(* ── electrode grouping & colours ── *)

splitByElectrode[goodSpecs_List] := Module[{groups},
  groups = GroupBy[goodSpecs, FileNameTake[#["File"], {-2}] &];
  {Keys[groups], Values[groups]}
]

electrodeStyles[n_Integer] :=
  Directive[ColorData[97][#], AbsoluteThickness[1.8]] & /@ Range[n]

electrodeLegend[electrodes_List, styles_List, th_Association] :=
  Placed[
    LineLegend[styles, electrodes,
      LegendMarkerSize -> 18,
      LabelStyle       -> Directive[th["Fg"], 12],
      Background       -> th["Bg"],
      LegendFunction   -> (Framed[#, Background -> th["Bg"],
        FrameStyle -> Directive[th["Fg"], AbsoluteThickness[1.1]]] &)
    ],
    Right
  ]

(* ── shared panel builder ── *)

extractTimeHr[specs_] := (#Spec["FinishTimeS"] / 3600. &) /@ specs

timeSeriesPanelMulti[
    tracesList_, styles_List, yLabel_, th_Association,
    logY_ : False, legend_ : None, opts___] :=
  Module[{plotFn, fg = th["Fg"], plotOpts, allPts, tLo, tHi, yHi},
    plotFn   = If[logY, ListLogPlot, ListLinePlot];
    plotOpts = FilterRules[{opts}, Options[plotFn]];
    allPts   = Flatten[tracesList, 1];
    {tLo, tHi} = If[allPts === {}, {0., 1.}, MinMax[allPts[[All, 1]]]];
    yHi        = If[allPts === {}, 1., Max[allPts[[All, 2]]] * 1.05];
    plotFn[
      tracesList,
      plotOpts,
      Joined      -> True,
      PlotMarkers -> {Automatic, 7},
      PlotStyle   -> styles,
      Sequence @@ ThemeChrome[th, 14, 1.1],
      FrameTicks  -> {{ThemeLinTicks[0., yHi, fg], None},
                      {ThemeLinTicks[tLo, tHi, fg], None}},
      FrameTicksStyle -> Directive[fg, 12],
      FrameLabel  -> {Style["Time (hours)", 14, fg], Style[yLabel, 14, fg]},
      PlotRange   -> {Automatic, {0, yHi}},
      PlotRangePadding -> Scaled[0.05],
      ImageSize   -> tsImageSize,
      AspectRatio -> tsAspect,
      Sequence @@ If[legend =!= None, {PlotLegends -> legend}, {}]
    ]
  ]

(* ── Rs ── *)

Options[PlotRsVsTime] = {Theme -> Automatic};

PlotRsVsTime[goodSpecs_List, opts : OptionsPattern[]] := Module[
  {th, electrodes, groups, n, styles, traces},
  th = ThemeFromOpts[{opts}];
  {electrodes, groups} = splitByElectrode[goodSpecs];
  n      = Length[electrodes];
  styles = electrodeStyles[n];
  traces = Table[
    Transpose[{extractTimeHr[groups[[i]]], (#Spec["Rs"] &) /@ groups[[i]]}],
    {i, n}];
  timeSeriesPanelMulti[traces, styles, "Rs (\[CapitalOmega])", th, False,
    electrodeLegend[electrodes, styles, th], opts]
]

(* ── C0 ── *)

Options[PlotC0VsTime] = {Theme -> Automatic};

PlotC0VsTime[goodSpecs_List, opts : OptionsPattern[]] := Module[
  {th, electrodes, groups, n, styles, traces},
  th = ThemeFromOpts[{opts}];
  {electrodes, groups} = splitByElectrode[goodSpecs];
  n      = Length[electrodes];
  styles = electrodeStyles[n];
  traces = Table[
    Transpose[{extractTimeHr[groups[[i]]], (1.*^9 #Spec["C0"] &) /@ groups[[i]]}],
    {i, n}];
  timeSeriesPanelMulti[traces, styles, "C0 (nF)", th, False,
    electrodeLegend[electrodes, styles, th], opts]
]

(* ── k_peak ── *)

(* Find peak k within a specific window; returns Missing if nothing found *)
peakInRange[spec_Association, kMin_?NumericQ, kMax_?NumericQ] := Module[
  {tau, g, k, mask, gM, kM},
  tau  = spec["Tau"];
  g    = spec["g"];
  k    = 1. / tau;
  mask = (kMin <= # <= kMax) & /@ k;
  gM   = Pick[g, mask];
  kM   = Pick[k, mask];
  If[Length[gM] < 1 || Max[gM] == 0.,
    Missing["NoPeak"],
    kM[[ First @ Ordering[gM, -1] ]]
  ]
]

(* per-electrode k_peak panel (linear y, outward ticks) *)
kPeakPanel[traces_, styles_, yLabel_, th_Association, yRange_, legend_, opts___] :=
  Module[{fg = th["Fg"], plotOpts, allK, allT, kLo, kHi, tLo, tHi},
    plotOpts = FilterRules[{opts}, Options[ListLinePlot]];
    allK = Flatten[traces[[All, All, 2]]];
    allT = Flatten[traces[[All, All, 1]]];
    {kLo, kHi} = If[yRange === All, If[allK === {}, {1., 10.}, MinMax[allK]], yRange];
    {tLo, tHi} = If[allT === {}, {0., 1.}, MinMax[allT]];
    ListLinePlot[
      traces,
      plotOpts,
      Joined      -> True,
      PlotMarkers -> {Automatic, 7},
      PlotStyle   -> styles,
      Sequence @@ ThemeChrome[th, 14, 1.1],
      FrameTicks  -> {{ThemeLinTicks[kLo, kHi, fg], None},
                      {ThemeLinTicks[tLo, tHi, fg], None}},
      FrameTicksStyle -> Directive[fg, 12],
      FrameLabel  -> {Style["Time (hours)", 14, fg], yLabel},
      PlotRange   -> {Automatic, {kLo, kHi}},
      PlotRangePadding -> Scaled[0.05],
      ImageSize   -> tsImageSize,
      AspectRatio -> tsAspect,
      Sequence @@ If[legend =!= None, {PlotLegends -> legend}, {}]
    ]
  ]

(* Range-constrained overload: one panel per k-range, per-electrode colouring *)
Options[PlotKPeakVsTime] = {Theme -> Automatic};

PlotKPeakVsTime[goodSpecs_List, kRanges_List, opts : OptionsPattern[]] := Module[
  {th, fg, electrodes, groups, n, styles, legend, panels, kMin, kMax, traces, tHr, kVals},
  th = ThemeFromOpts[{opts}];
  fg = th["Fg"];
  {electrodes, groups} = splitByElectrode[goodSpecs];
  n      = Length[electrodes];
  styles = electrodeStyles[n];
  legend = electrodeLegend[electrodes, styles, th];

  panels = Table[
    kMin = kRanges[[r, 1]];  kMax = kRanges[[r, 2]];
    traces = Table[
      tHr   = extractTimeHr[groups[[i]]];
      kVals = peakInRange[#, kMin, kMax] & /@ groups[[i, All, "Spec"]];
      Select[Transpose[{tHr, kVals}], NumericQ[#[[2]]] &],
      {i, n}];
    kPeakPanel[traces, styles,
      Style[Row[{"k", Subscript["peak", ""], "  [",
        kMin, "\[Dash]", kMax, "] (", Superscript["s", -1], ")"}], 13, fg],
      th, {kMin, kMax}, If[r == 1, legend, None], opts],
    {r, Length[kRanges]}
  ];

  Row[panels, Spacer[20]]
]

PlotKPeakVsTime[goodSpecs_List, opts : OptionsPattern[]] := Module[
  {th, fg, electrodes, groups, n, styles, traces, tHr, peaks, kVals},
  th = ThemeFromOpts[{opts}];
  fg = th["Fg"];
  {electrodes, groups} = splitByElectrode[goodSpecs];
  n      = Length[electrodes];
  styles = electrodeStyles[n];
  traces = Table[
    tHr   = extractTimeHr[groups[[i]]];
    peaks = DCTPeakAnalysis`FindDCTPeak /@ (groups[[i, All, "Spec"]]);
    kVals = Map[If[AssociationQ[#], #["kPeak"], Missing["NoPeak"]] &, peaks];
    Select[Transpose[{tHr, kVals}], NumericQ[#[[2]]] &],
    {i, n}];
  kPeakPanel[traces, styles,
    Style[Row[{"k", Subscript["", "peak"], " (", Superscript["s", -1], ")"}], 14, fg],
    th, All, electrodeLegend[electrodes, styles, th], opts]
]

(* ── combined grid (legend on Rs panel only) ── *)

Options[PlotTimeSeries] = {Theme -> Automatic};

PlotTimeSeries[goodSpecs_List, opts : OptionsPattern[]] := Module[
  {th, fg, electrodes, groups, n, styles, legend,
   rsTraces, c0Traces, kTraces, tHr, peaks, kVals,
   rsPlot, c0Plot, kPlot},

  th = ThemeFromOpts[{opts}];
  fg = th["Fg"];
  {electrodes, groups} = splitByElectrode[goodSpecs];
  n      = Length[electrodes];
  styles = electrodeStyles[n];
  legend = electrodeLegend[electrodes, styles, th];

  rsTraces = Table[
    Transpose[{extractTimeHr[groups[[i]]], (#Spec["Rs"] &) /@ groups[[i]]}],
    {i, n}];
  c0Traces = Table[
    Transpose[{extractTimeHr[groups[[i]]], (#Spec["C0"] &) /@ groups[[i]]}],
    {i, n}];
  kTraces = Table[
    tHr   = extractTimeHr[groups[[i]]];
    peaks = DCTPeakAnalysis`FindDCTPeak /@ (groups[[i, All, "Spec"]]);
    kVals = Map[If[AssociationQ[#], #["kPeak"], Missing["NoPeak"]] &, peaks];
    Select[Transpose[{tHr, kVals}], NumericQ[#[[2]]] &],
    {i, n}];

  (* panels without per-panel legends; one shared legend for the whole figure *)
  rsPlot = timeSeriesPanelMulti[rsTraces, styles, "Rs (\[CapitalOmega])", th, False, None];
  c0Plot = timeSeriesPanelMulti[c0Traces, styles, "C0 (F)", th, False, None];
  kPlot  = kPeakPanel[kTraces, styles,
    Style[Row[{"k", Subscript["", "peak"], " (", Superscript["s", -1], ")"}], 14, fg],
    th, All, None];

  (* One shared legend outside the grid (a per-panel legend shrank the Rs cell
     and clipped its y-label). *)
  Legended[
    GraphicsGrid[
      {{rsPlot, c0Plot, kPlot}},
      ImageSize  -> 1200,
      Spacings   -> {0.3, 0},
      Background -> th["Bg"]
    ],
    legend
  ]
]

End[]
EndPackage[]
