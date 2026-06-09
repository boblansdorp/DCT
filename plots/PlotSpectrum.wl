(* ::Package:: *)

(* DCTPlots`PlotSpectrum
   g(tau) and g(k) distribution plots for one or more DCT spectra.
   All functions are pure: take result Associations, return Graphics.

   Public API:
     PlotGK[goodSpecs]         rainbow overlay of g(k) for all spectra
     PlotGTau[goodSpecs]       rainbow overlay of g(tau)
     PlotCumulative[goodSpecs] cumulative capacitance vs k

   Default colour scheme is dark (white-on-dark) for wolfbook display.
   For light-background export, pass: Background->White,
     FrameStyle->Directive[Black,AbsoluteThickness[1.2]],
     LabelStyle->Directive[Black,16,FontFamily->"Arial"]
*)

BeginPackage["DCTPlots`"]

PlotGK::usage =
  "PlotGK[goodSpecs] plots g(k) (electron-transfer rate space) for all \
spectra in goodSpecs, rainbow-coloured. Optional second arg is a list of \
labels."

PlotGTau::usage =
  "PlotGTau[goodSpecs] plots g(tau) for all spectra, rainbow-coloured."

PlotCumulative::usage =
  "PlotCumulative[goodSpecs] plots cumulative capacitance int g dk vs k."

Begin["`Private`"]

rainbowColors[n_Integer] :=
  If[n <= 1,
    {ColorData["Rainbow"][0.2]},
    Table[ColorData["Rainbow"][u], {u, 0.05, 0.95, 0.90 / (n - 1)}]
  ]

defaultLabels[goodSpecs_List] :=
  FileNameTake /@ goodSpecs[[All, "File"]]

lineStyles[colors_List] :=
  Directive[#, AbsoluteThickness[2.2]] & /@ colors

$darkBg    = GrayLevel[0.12];
$fgStyle   = Directive[White, AbsoluteThickness[1.2]];
$lblStyle  = Directive[White, 16, FontFamily -> "Arial"];

(* ------------------------------------------------------------------ *)

PlotGK[goodSpecs_List, labels_ : Automatic, opts : OptionsPattern[]] := Module[
  {n, traces, cols, styles, lbls, kMin, kMax, xTicks, gMax},

  n      = Length[goodSpecs];
  traces = (Transpose[{1. / #Spec["Tau"], #Spec["g"]}] &) /@ goodSpecs;
  cols   = rainbowColors[n];
  styles = lineStyles[cols];
  lbls   = If[labels === Automatic, defaultLabels[goodSpecs], labels];

  kMin = Min[Flatten[traces[[All, All, 1]]]];
  kMax = Max[Flatten[traces[[All, All, 1]]]];
  xTicks = Table[
    {10.^p, Superscript[10, p]},
    {p, Floor[Log10[kMin]], Ceiling[Log10[kMax]]}
  ];

  (* 99th-percentile y cap: edge bins and outlier fits can spike far above the
     main spectral features and collapse everything else to a flat line. *)
  gMax = Quantile[Flatten[traces[[All, All, 2]]], 0.99] * 1.1;

  ListLinePlot[
    traces,
    opts,
    PlotStyle        -> styles,
    ScalingFunctions -> {"Log10", None},
    Frame            -> True,
    Axes             -> False,
    Background       -> $darkBg,
    FrameStyle       -> $fgStyle,
    FrameTicks       -> {{Automatic, None}, {xTicks, None}},
    FrameTicksStyle  -> Directive[White, 16],
    LabelStyle       -> $lblStyle,
    FrameLabel       -> {
      Style[Row[{"Electron-transfer rate, k (", Superscript["s", -1], ")"}], 16, White],
      Style["g(k) (F/decade)", 16, White]
    },
    PlotRange        -> {Automatic, {0, gMax}},
    PlotRangePadding -> {{Scaled[0.02], Scaled[0.02]}, {Scaled[0.02], Scaled[0.05]}},
    ImageSize        -> 600,
    AspectRatio      -> 0.65,
    PlotLegends      -> Placed[LineLegend[styles, lbls,
      LegendMarkerSize -> 24,
      LabelStyle       -> Directive[White, 13],
      Background       -> $darkBg,
      LegendFunction   -> (Framed[#, Background -> $darkBg, FrameStyle -> $fgStyle] &)],
      Right]
  ]
]


(* ------------------------------------------------------------------ *)

cumulativeVsK[spec_Association] := Module[
  {tau, g, k, logK, n, areas, cumC},
  tau  = spec["Tau"];
  g    = spec["g"];
  k    = 1. / tau;
  logK = Log10[k];
  n    = Length[g];
  If[n < 2, Return[{}]];
  (* trapezoid areas, always positive regardless of k sort order *)
  areas = Table[0.5 (g[[j]] + g[[j+1]]) Abs[logK[[j+1]] - logK[[j]]], {j, n-1}];
  (* integrate from k_min (left) upward: zero at low k, rises toward high k *)
  cumC = Append[Reverse @ Accumulate[Reverse[areas]], 0.];
  Transpose[{k, cumC}]
]

PlotCumulative[goodSpecs_List, labels_ : Automatic, opts : OptionsPattern[]] := Module[
  {n, traces, cols, styles, lbls},

  n      = Length[goodSpecs];
  traces = cumulativeVsK /@ (goodSpecs[[All, "Spec"]]);
  cols   = rainbowColors[n];
  styles = lineStyles[cols];
  lbls   = If[labels === Automatic, defaultLabels[goodSpecs], labels];

  ListLinePlot[
    traces,
    PlotStyle        -> styles,
    ScalingFunctions -> {"Log10", None},
    Frame            -> True,
    Axes             -> False,
    Background       -> $darkBg,
    FrameStyle       -> $fgStyle,
    LabelStyle       -> $lblStyle,
    FrameLabel       -> {
      Style[Row[{"k (", Superscript["s", -1], ")"}], 16, White],
      Style["Cumulative capacitance (F)", 16, White]
    },
    PlotRange        -> All,
    PlotRangePadding -> Scaled[0.04],
    ImageSize        -> 600,
    AspectRatio      -> 0.65,
    PlotLegends      -> Placed[LineLegend[styles, lbls,
      LegendMarkerSize -> 24,
      LabelStyle       -> Directive[White, 13],
      Background       -> $darkBg,
      LegendFunction   -> (Framed[#, Background -> $darkBg, FrameStyle -> $fgStyle] &)],
      Right],
    opts
  ]
]

End[]
EndPackage[]
