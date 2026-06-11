(* ::Package:: *)

(* DCTPlots`PlotSpectrum
   g(tau) and g(k) distribution plots for one or more DCT spectra.
   All functions are pure: take result Associations, return Graphics.

   Public API:
     PlotGK[goodSpecs]         rainbow overlay of g(k) for all spectra
     PlotGTau[goodSpecs]       rainbow overlay of g(tau)
     PlotCumulative[goodSpecs] cumulative capacitance vs k

   Log x-axis is drawn with native ListLogLinearPlot (Mathematica picks the
   range/scaling); the theme sets colours + outward FrameTicks (ThemeLogTicks /
   ThemeLinTicks). Default colour scheme is dark (white-on-dark); pass
   Theme -> "Publication" for a white-background figure.
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

HourMarkLabels::usage =
  "HourMarkLabels[specs] / HourMarkLabels[specs, step] returns a legend-label \
list: the spectrum nearest each 0, step, 2 step, ... hour mark (step default 6) \
is labelled \"0 h\", \"6 h\", ...; every other entry is \"\" (blank). Pair with \
PlotGK so the drift legend shows clean round-hour marks, not file names."

Begin["`Private`"]

rainbowColors[n_Integer] :=
  If[n <= 1,
    {ColorData["Rainbow"][0.2]},
    Table[ColorData["Rainbow"][u], {u, 0.05, 0.95, 0.90 / (n - 1)}]
  ]

defaultLabels[goodSpecs_List] :=
  FileNameTake /@ goodSpecs[[All, "File"]]

HourMarkLabels[specs_List, step_ : 6] := Module[
  {t0, hrs, marks, idx, labels},
  t0     = Min[(#["Spec"]["FinishTimeS"]) & /@ specs];
  hrs    = ((#["Spec"]["FinishTimeS"] - t0) / 3600.) & /@ specs;
  marks  = Range[0, step * Round[Max[hrs] / step], step];
  idx    = (First @ Ordering[Abs[hrs - #], 1]) & /@ marks;
  labels = ConstantArray["", Length[specs]];
  Do[labels[[idx[[i]]]] = ToString[marks[[i]]] <> " h", {i, Length[marks]}];
  labels
]

lineStyles[colors_List] :=
  Directive[#, AbsoluteThickness[2.2]] & /@ colors

(* Build the legend entry list: keep only curves with a non-blank label (callers
   can blank out the rest, e.g. HourMarkLabels), then cap at maxN evenly-spaced
   so a per-file legend can't get enormous. All traces stay plotted regardless. *)
legendSubset[styles_List, lbls_List, maxN_Integer : 8] := Module[{keep, s, l},
  keep = Select[Range[Length[lbls]], (lbls[[#]] =!= "" && lbls[[#]] =!= Null) &];
  If[keep === {}, keep = Range[Length[lbls]]];
  {s, l} = {styles[[keep]], lbls[[keep]]};
  If[Length[s] <= maxN,
    {s, l},
    With[{idx = DeleteDuplicates @ Round @ Subdivide[1, Length[s], maxN - 1]},
      {s[[idx]], l[[idx]]}]
  ]
]

themedLegend[styles_, lbls_, th_, fg_] :=
  Placed[LineLegend[styles, lbls,
    LegendMarkerSize -> 24,
    LabelStyle       -> Directive[fg, 13],
    Background       -> th["Bg"],
    LegendFunction   -> (Framed[#, Background -> th["Bg"],
      FrameStyle -> Directive[fg, AbsoluteThickness[1.2]]] &)],
    Right]

(* ------------------------------------------------------------------ *)

Options[PlotGK] = {Theme -> Automatic, "GScale" -> 1, "AreaNorm" -> False};

PlotGK[goodSpecs_List, labels : (_List | Automatic) : Automatic,
       opts : OptionsPattern[]] := Module[
  {th, fg, plotOpts, gScale, areaNorm, gLabel, n, traces, cols, styles, lbls,
   legStyles, legLbls, kMin, kMax, gPeak, gMax},

  th       = ThemeFromOpts[{opts}];   (* no OptionValue -> no nodef on pass-throughs *)
  fg       = th["Fg"];
  plotOpts = FilterRules[{opts}, Options[ListLogLinearPlot]];  (* drop Theme; user opts win *)

  (* string-keyed options (a bare-symbol option would resolve to Global` in the
     notebook). GScale rescales g (e.g. 10^9 -> nF); AreaNorm divides each curve
     by its integrated capacitance. *)
  gScale   = "GScale" /. {opts} /. "GScale" -> 1;
  areaNorm = TrueQ["AreaNorm" /. {opts} /. "AreaNorm" -> False];

  n      = Length[goodSpecs];
  traces = (If[areaNorm,
      Transpose[{1. / #Spec["Tau"],
        #Spec["g"] / Max[DCTPeakAnalysis`IntegrateSpectrum[#Spec], 1.*^-15]}],
      Transpose[{1. / #Spec["Tau"], gScale #Spec["g"]}]
    ] &) /@ goodSpecs;
  cols   = rainbowColors[n];
  styles = lineStyles[cols];
  lbls   = If[labels === Automatic, defaultLabels[goodSpecs], labels];
  {legStyles, legLbls} = legendSubset[styles, lbls];

  kMin = Min[Flatten[traces[[All, All, 1]]]];
  kMax = Max[Flatten[traces[[All, All, 1]]]];

  (* y cap from the 98th percentile of per-spectrum PEAK heights (not of all g
     values) so an outlier fit can't blow up the axis. *)
  gPeak = Quantile[Max /@ traces[[All, All, 2]], 0.98];
  gMax  = gPeak * 1.1;

  gLabel = Which[
    areaNorm,        Row[{"g(k) / C", Subscript["dist", ""], "  (", Superscript["decade", -1], ")"}],
    gScale == 10.^9, "g(k) (nF/decade)",
    True,            "g(k) (F/decade)"
  ];

  ListLogLinearPlot[
    traces,
    plotOpts,                                  (* user overrides win (first) *)
    Joined      -> True,
    PlotStyle   -> styles,
    Sequence @@ ThemeChrome[th, 16, 1.2],
    FrameTicks       -> {{ThemeLinTicks[0., gMax, fg], None},
                         {ThemeLogTicks[kMin, kMax, fg], None}},
    FrameTicksStyle  -> Directive[fg, 16],
    FrameLabel  -> {
      Style[Row[{"Electron-transfer rate, k (", Superscript["s", -1], ")"}], 16, fg],
      Style[gLabel, 16, fg]
    },
    PlotRange        -> {All, {0, gMax}},
    PlotRangePadding -> {{Scaled[0.02], Scaled[0.02]}, {Scaled[0.02], 0}},
    ImageSize        -> 600,
    AspectRatio      -> 0.65,
    PlotLegends      -> themedLegend[legStyles, legLbls, th, fg]
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

Options[PlotCumulative] = {Theme -> Automatic};

PlotCumulative[goodSpecs_List, labels : (_List | Automatic) : Automatic,
               opts : OptionsPattern[]] := Module[
  {th, fg, plotOpts, n, traces, cols, styles, lbls, legStyles, legLbls,
   kMin, kMax, cMax},

  th       = ThemeFromOpts[{opts}];
  fg       = th["Fg"];
  plotOpts = FilterRules[{opts}, Options[ListLogLinearPlot]];

  n      = Length[goodSpecs];
  traces = cumulativeVsK /@ (goodSpecs[[All, "Spec"]]);
  cols   = rainbowColors[n];
  styles = lineStyles[cols];
  lbls   = If[labels === Automatic, defaultLabels[goodSpecs], labels];
  {legStyles, legLbls} = legendSubset[styles, lbls];

  kMin = Min[Select[Flatten[traces[[All, All, 1]]], Positive]];
  kMax = Max[Flatten[traces[[All, All, 1]]]];
  cMax = Max[Flatten[traces[[All, All, 2]]]] * 1.05;

  ListLogLinearPlot[
    traces,
    plotOpts,
    Joined      -> True,
    PlotStyle   -> styles,
    Sequence @@ ThemeChrome[th, 16, 1.2],
    FrameTicks       -> {{ThemeLinTicks[0., cMax, fg], None},
                         {ThemeLogTicks[kMin, kMax, fg], None}},
    FrameTicksStyle  -> Directive[fg, 16],
    FrameLabel  -> {
      Style[Row[{"k (", Superscript["s", -1], ")"}], 16, fg],
      Style["Cumulative capacitance (F)", 16, fg]
    },
    PlotRange        -> {All, {0, cMax}},
    PlotRangePadding -> Scaled[0.04],
    ImageSize        -> 600,
    AspectRatio      -> 0.65,
    PlotLegends      -> themedLegend[legStyles, legLbls, th, fg]
  ]
]

End[]
EndPackage[]
