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

   Default colour scheme is dark (white-on-dark) for wolfbook display.
   For light-background export pass: Background->White,
     FrameStyle->Directive[Black,AbsoluteThickness[1.1]],
     LabelStyle->Directive[Black,14,FontFamily->"Arial"]
*)

BeginPackage["DCTPlots`", {"DCTPeakAnalysis`"}]

PlotRsVsTime::usage =
  "PlotRsVsTime[goodSpecs] plots Rs (\[CapitalOmega]) vs time (hours), \
one trace per electrode."

PlotC0VsTime::usage =
  "PlotC0VsTime[goodSpecs] plots C0 (F) vs time (hours), \
one trace per electrode."

PlotKPeakVsTime::usage =
  "PlotKPeakVsTime[goodSpecs] plots the global dominant peak k (s^-1) vs time. \
PlotKPeakVsTime[goodSpecs, {{kMin1,kMax1},{kMin2,kMax2},...}] shows one panel \
per k-range, tracking the peak within each window."

PlotTimeSeries::usage =
  "PlotTimeSeries[goodSpecs] returns a GraphicsGrid with Rs, C0, and kPeak \
vs time panels, each electrode in a distinct colour."

Begin["`Private`"]

$darkBg      = GrayLevel[0.12];
tsFrameStyle = Directive[White, AbsoluteThickness[1.1]];
tsLabelStyle = Directive[White, 14, FontFamily -> "Arial"];
tsImageSize  = 380;
tsAspect     = 0.75;

(* ── electrode grouping & colours ── *)

splitByElectrode[goodSpecs_List] := Module[{groups},
  groups = GroupBy[goodSpecs, FileNameTake[#["File"], {-2}] &];
  {Keys[groups], Values[groups]}
]

electrodeStyles[n_Integer] :=
  Directive[ColorData[97][#], AbsoluteThickness[1.8]] & /@ Range[n]

electrodeLegend[electrodes_List, styles_List] :=
  Placed[
    LineLegend[styles, electrodes,
      LegendMarkerSize -> 18,
      LabelStyle       -> Directive[White, 12],
      Background       -> $darkBg,
      LegendFunction   -> (Framed[#, Background -> $darkBg,
                              FrameStyle -> tsFrameStyle] &)
    ],
    Right
  ]

(* ── shared panel builder ── *)

extractTimeHr[specs_] := (#Spec["FinishTimeS"] / 3600. &) /@ specs

timeSeriesPanelMulti[
    tracesList_, styles_List, yLabel_,
    logY_ : False, legend_ : None, opts___] :=
  Module[{plotFn},
    plotFn = If[logY, ListLogPlot, ListLinePlot];
    plotFn[
      tracesList,
      Joined      -> True,
      PlotMarkers -> {Automatic, 7},
      PlotStyle   -> styles,
      Frame       -> True,
      Axes        -> False,
      Background  -> $darkBg,
      FrameStyle  -> tsFrameStyle,
      LabelStyle  -> tsLabelStyle,
      FrameLabel  -> {Style["Time (hours)", 14, White], Style[yLabel, 14, White]},
      PlotRange   -> {Automatic, {0, Automatic}},
      PlotRangePadding -> Scaled[0.05],
      ImageSize   -> tsImageSize,
      AspectRatio -> tsAspect,
      Sequence @@ If[legend =!= None, {PlotLegends -> legend}, {}],
      opts
    ]
  ]

(* ── Rs ── *)

PlotRsVsTime[goodSpecs_List, opts : OptionsPattern[]] := Module[
  {electrodes, groups, n, styles, traces},
  {electrodes, groups} = splitByElectrode[goodSpecs];
  n      = Length[electrodes];
  styles = electrodeStyles[n];
  traces = Table[
    Transpose[{extractTimeHr[groups[[i]]], (#Spec["Rs"] &) /@ groups[[i]]}],
    {i, n}];
  timeSeriesPanelMulti[traces, styles, "Rs (\[CapitalOmega])", False,
    electrodeLegend[electrodes, styles], opts]
]

(* ── C0 ── *)

PlotC0VsTime[goodSpecs_List, opts : OptionsPattern[]] := Module[
  {electrodes, groups, n, styles, traces},
  {electrodes, groups} = splitByElectrode[goodSpecs];
  n      = Length[electrodes];
  styles = electrodeStyles[n];
  traces = Table[
    Transpose[{extractTimeHr[groups[[i]]], (#Spec["C0"] &) /@ groups[[i]]}],
    {i, n}];
  timeSeriesPanelMulti[traces, styles, "C0 (F)", False,
    electrodeLegend[electrodes, styles], opts]
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

(* Range-constrained overload: one panel per k-range, per-electrode colouring *)
PlotKPeakVsTime[goodSpecs_List, kRanges_List, opts : OptionsPattern[]] := Module[
  {electrodes, groups, n, styles, legend, panels, kMin, kMax, traces, tHr, kVals},
  {electrodes, groups} = splitByElectrode[goodSpecs];
  n      = Length[electrodes];
  styles = electrodeStyles[n];
  legend = electrodeLegend[electrodes, styles];

  panels = Table[
    kMin = kRanges[[r, 1]];  kMax = kRanges[[r, 2]];
    traces = Table[
      tHr   = extractTimeHr[groups[[i]]];
      kVals = peakInRange[#, kMin, kMax] & /@ groups[[i, All, "Spec"]];
      Select[Transpose[{tHr, kVals}], NumericQ[#[[2]]] &],
      {i, n}];
    ListLogPlot[
      traces,
      Joined      -> True,
      PlotMarkers -> {Automatic, 7},
      PlotStyle   -> styles,
      Frame       -> True,
      Axes        -> False,
      Background  -> $darkBg,
      FrameStyle  -> tsFrameStyle,
      LabelStyle  -> tsLabelStyle,
      FrameLabel  -> {
        Style["Time (hours)", 14, White],
        Style[Row[{"k", Subscript["peak", ""], "  [",
          kMin, "\[Dash]", kMax, "] (", Superscript["s", -1], ")"}], 13, White]
      },
      PlotRange        -> {Automatic, {kMin, kMax}},
      PlotRangePadding -> Scaled[0.05],
      ImageSize        -> tsImageSize,
      AspectRatio      -> tsAspect,
      Sequence @@ If[r == 1, {PlotLegends -> legend}, {}],
      opts
    ],
    {r, Length[kRanges]}
  ];

  Row[panels, Spacer[20]]
]

PlotKPeakVsTime[goodSpecs_List, opts : OptionsPattern[]] := Module[
  {electrodes, groups, n, styles, traces, tHr, peaks, kVals},
  {electrodes, groups} = splitByElectrode[goodSpecs];
  n      = Length[electrodes];
  styles = electrodeStyles[n];
  traces = Table[
    tHr   = extractTimeHr[groups[[i]]];
    peaks = DCTPeakAnalysis`FindDCTPeak /@ (groups[[i, All, "Spec"]]);
    kVals = Map[If[AssociationQ[#], #["kPeak"], Missing["NoPeak"]] &, peaks];
    Select[Transpose[{tHr, kVals}], NumericQ[#[[2]]] &],
    {i, n}];
  ListLogPlot[
    traces,
    Joined      -> True,
    PlotMarkers -> {Automatic, 7},
    PlotStyle   -> styles,
    Frame       -> True,
    Axes        -> False,
    Background  -> $darkBg,
    FrameStyle  -> tsFrameStyle,
    LabelStyle  -> tsLabelStyle,
    FrameLabel  -> {
      Style["Time (hours)", 14, White],
      Style[Row[{"k", Subscript["", "peak"], " (", Superscript["s", -1], ")"}], 14, White]
    },
    PlotRange        -> All,
    PlotRangePadding -> Scaled[0.05],
    ImageSize        -> tsImageSize,
    AspectRatio      -> tsAspect,
    PlotLegends      -> electrodeLegend[electrodes, styles],
    opts
  ]
]

(* ── combined grid (legend on Rs panel only) ── *)

PlotTimeSeries[goodSpecs_List, opts : OptionsPattern[]] := Module[
  {electrodes, groups, n, styles, legend,
   rsTraces, c0Traces, kTraces, tHr, peaks, kVals,
   rsPlot, c0Plot, kPlot},

  {electrodes, groups} = splitByElectrode[goodSpecs];
  n      = Length[electrodes];
  styles = electrodeStyles[n];
  legend = electrodeLegend[electrodes, styles];

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

  (* legend on Rs panel only to avoid repetition *)
  rsPlot = timeSeriesPanelMulti[rsTraces, styles, "Rs (\[CapitalOmega])", False, legend];
  c0Plot = timeSeriesPanelMulti[c0Traces, styles, "C0 (F)", False, None];
  kPlot  = ListLogPlot[
    kTraces,
    Joined      -> True,
    PlotMarkers -> {Automatic, 7},
    PlotStyle   -> styles,
    Frame       -> True,
    Axes        -> False,
    Background  -> $darkBg,
    FrameStyle  -> tsFrameStyle,
    LabelStyle  -> tsLabelStyle,
    FrameLabel  -> {
      Style["Time (hours)", 14, White],
      Style[Row[{"k", Subscript["", "peak"], " (", Superscript["s", -1], ")"}], 14, White]
    },
    PlotRange        -> All,
    PlotRangePadding -> Scaled[0.05],
    ImageSize        -> tsImageSize,
    AspectRatio      -> tsAspect];

  GraphicsGrid[
    {{rsPlot, c0Plot, kPlot}},
    ImageSize  -> 1200,
    Spacings   -> {0.4, 0.4},
    Background -> $darkBg,
    opts
  ]
]

End[]
EndPackage[]
