(* ::Package:: *)

(* DCTPlots — admittance statistics visualisation
   Two plots for Section 2 of the analysis notebook.

   Public API:
     PlotCV[stats]
     PlotCV[stats, cvThresh]
       Log-x plot of CV(f) = sigma/mean vs frequency.
       Red dashed line marks cvThresh (default 0.1).

     PlotVariancePowerLaw[stats]
       Log-log scatter of Var(|Y(f)|) vs |Y_mean(f)|.
       Points coloured by log10(f), best-fit power law overlaid.
       Derived alpha and weightPower printed in the plot.

   Theme-aware: default DCTPlots`$DCTTheme ("Dark" | "Publication"); override
   one call with Theme -> "Publication". Rainbow point colours and the red
   threshold line are theme-independent.
*)

BeginPackage["DCTPlots`"]

PlotCV::usage =
  "PlotCV[stats] or PlotCV[stats, cvThresh] plots CV(f) = sigma/mean of \
|Y(f)| vs frequency (log x). Red dashed line at cvThresh (default 0.1)."

PlotVariancePowerLaw::usage =
  "PlotVariancePowerLaw[stats] or PlotVariancePowerLaw[stats, fMin, fMax] \
plots Var(|Y(f)|) vs |Y_mean(f)| on log-log axes, coloured by log10(f), \
with the best-fit power law overlaid. The 3-arg form fits only in [fMin,fMax] \
and marks the boundaries with dashed vertical lines; the scatter shows all data."

Begin["`Private`"]

(* ------------------------------------------------------------------ *)

Options[PlotCV] = {Theme -> Automatic};

PlotCV[stats_Association, cvThresh : (_?NumericQ) : 0.1,
       opts : OptionsPattern[]] := Module[
  {th, fg, plotOpts, freq, cv, cvPairs, threshPair, yMax, fMin, fMax,
   cvLineStyle, threshStyle},

  th       = ThemeFromOpts[{opts}];
  fg       = th["Fg"];
  plotOpts = FilterRules[{opts}, Options[ListLinePlot]];

  freq = stats["FreqHz"];
  cv   = stats["CV"];

  cvPairs      = Transpose[{freq, cv}];
  threshPair   = {{Min[freq], cvThresh}, {Max[freq], cvThresh}};
  {fMin, fMax} = MinMax[freq];
  yMax         = Max[Select[cv, NumericQ]] * 1.1;
  cvLineStyle  = Directive[ColorData["Rainbow"][0.55], AbsoluteThickness[2.2]];
  threshStyle  = Directive[Red, Dashed, AbsoluteThickness[1.8]];

  ListLinePlot[
    {cvPairs, threshPair},
    plotOpts,
    ScalingFunctions  -> {"Log10", None},
    PlotStyle         -> {cvLineStyle, threshStyle},
    Sequence @@ ThemeChrome[th, 14, 1.2],
    FrameTicks        -> {{ThemeLinTicks[0., yMax, fg], None},
                          {ThemeLogTicks[fMin, fMax, fg], None}},
    FrameTicksStyle   -> Directive[fg, 13],
    FrameLabel        -> {
      Style["Frequency (Hz)", 14, fg],
      Style["\[Sigma]/\[Mu] of |Y(f)|", 14, fg]
    },
    PlotRange         -> {{fMin, fMax}, {0, yMax}},
    PlotRangePadding  -> {{Scaled[0.02], Scaled[0.02]}, {0, Scaled[0.05]}},
    ImageSize         -> 600,
    AspectRatio       -> 0.55,
    PlotLegends       -> Placed[
      LineLegend[
        {cvLineStyle, threshStyle},
        {"CV(f)", "Threshold " <> ToString[cvThresh]},
        LabelStyle      -> Directive[fg, 12],
        Background      -> th["Bg"],
        LegendFunction  -> (Framed[#, Background -> th["Bg"],
                              FrameStyle -> Directive[fg, AbsoluteThickness[1.2]]] &)
      ],
      {Right, Top}
    ]
  ]
]

(* ------------------------------------------------------------------ *)

(* Shared Show builder for the variance scatter + fit + epilog primitives.
   Built on ScalingFunctions -> {"Log10","Log10"} so the theme outward-tick
   generators apply; the scaled coordinate system is log10, so epilog
   primitives are placed in Log10 (not natural-log) coordinates. xRange/yRange
   are the real-value data extents used for the outward log ticks. *)
variancePowerLawShow[th_Association, singletonSets_, singletonStyles_,
    fitLine_, xRange_, yRange_, epilogPrims_] :=
  Module[{fg = th["Fg"]},
    Show[
      ListPlot[
        singletonSets,
        PlotTheme         -> "Default",
        ScalingFunctions  -> {"Log10", "Log10"},
        PlotStyle         -> singletonStyles,
        Frame             -> True,
        Axes              -> False,
        Background        -> th["Bg"],
        FrameStyle        -> Directive[fg, AbsoluteThickness[1.2]],
        FrameTicks        -> {{ThemeLogTicks[yRange[[1]], yRange[[2]], fg], None},
                              {ThemeLogTicks[xRange[[1]], xRange[[2]], fg], None}},
        FrameTicksStyle   -> Directive[fg, 13],
        LabelStyle        -> Directive[fg, 14, FontFamily -> "Arial"],
        FrameLabel        -> {
          {Style["Var(|Y(f)|) (S^2)", 14, fg], None},
          {Style["\[LeftAngleBracket]|Y(f)|\[RightAngleBracket] (S)", 14, fg], None}
        },
        PlotRangePadding  -> Scaled[0.04],
        ImageSize         -> 600,
        AspectRatio       -> 0.65
      ],
      ListLinePlot[
        {fitLine},
        PlotTheme        -> "Default",
        ScalingFunctions -> {"Log10", "Log10"},
        PlotStyle        -> Directive[fg, AbsoluteThickness[2.2]]
      ],
      Epilog -> epilogPrims
    ]
  ]

annotLabelFor[alpha_] := Row[{
  "Var \[Proportional] |Y|",
  Superscript["\[Alpha]", ""],
  ",  \[Alpha] = ", NumberForm[alpha, {4, 2}],
  "   \[DoubleRightArrow]   weightPower = ",
  NumberForm[-alpha, {4, 2}]
}]

(* frequency decade labels near data points (Log10 = scaled-axis coords) *)
freqLabelsFor[freqGood_, yMeanGood_, yVarGood_, fg_] := Module[
  {freqTargets = {1., 10., 100., 1000., 10000.}, pts},
  pts = Map[
    Function[f,
      Module[{idx, logDiff},
        idx     = First @ Ordering[Abs[Log10[freqGood] - Log10[f]], 1];
        logDiff = Abs[Log10[freqGood[[idx]]] - Log10[f]];
        If[logDiff < 0.3, {yMeanGood[[idx]], yVarGood[[idx]], f}, Nothing]
      ]
    ],
    freqTargets
  ];
  Style[Text[Row[{ToString[Round[#[[3]]]], " Hz"}],
    {Log10[#[[1]]], Log10[#[[2]]]}, {-1.3, 0.6}], fg, 10] & /@ pts
]

(* ------------------------------------------------------------------ *)

Options[PlotVariancePowerLaw] = {Theme -> Automatic};

PlotVariancePowerLaw[stats_Association, opts : OptionsPattern[]] := Module[
  {th, fg, freq, yMean, yVar, alpha, c, logFreq, logFMin, logFMax,
   goodMask, yMeanGood, yVarGood, freqGood,
   colsGood, singletonSets, singletonStyles,
   xFitVals, yFitVals, fitLine, xRange, yRange, epilog},

  th = ThemeFromOpts[{opts}];
  fg = th["Fg"];

  freq  = stats["FreqHz"];
  yMean = stats["YMean"];
  yVar  = stats["YVar"];
  alpha = stats["PowerLawAlpha"];
  c     = stats["PowerLawC"];

  logFreq = Log10[freq];
  logFMin = Min[logFreq]; logFMax = Max[logFreq];

  goodMask  = MapThread[#1 > 0 && #2 > 0 &, {yMean, yVar}];
  yMeanGood = Pick[yMean, goodMask];
  yVarGood  = Pick[yVar,  goodMask];
  freqGood  = Pick[freq,  goodMask];

  colsGood = ColorData["Rainbow"] /@ Rescale[Log10[freqGood], {logFMin, logFMax}];
  singletonSets   = {{#1, #2}} & @@@ Transpose[{yMeanGood, yVarGood}];
  singletonStyles = Directive[#, PointSize[0.013]] & /@ colsGood;

  xFitVals = 10.^Range[Log10[Min[yMeanGood]] - 0.05, Log10[Max[yMeanGood]] + 0.05, 0.02];
  yFitVals = c * xFitVals^alpha;
  fitLine  = Transpose[{xFitVals, yFitVals}];

  xRange = MinMax[yMeanGood];
  yRange = MinMax[yVarGood];

  epilog = Join[
    {Style[Text[annotLabelFor[alpha], Scaled[{0.04, 0.96}], {-1, 1}], fg, 12]},
    freqLabelsFor[freqGood, yMeanGood, yVarGood, fg]
  ];

  variancePowerLawShow[th, singletonSets, singletonStyles, fitLine, xRange, yRange, epilog]
]

(* ------------------------------------------------------------------ *)

PlotVariancePowerLaw[stats_Association, fMin_?NumericQ, fMax_?NumericQ,
                     opts : OptionsPattern[]] := Module[
  {th, fg, freq, yMean, yVar, logFreq, logFMin, logFMax,
   goodMask, yMeanGood, yVarGood, freqGood,
   fitMask, yMeanFit, yVarFit, colsGood, singletonSets, singletonStyles,
   goodPairs, goodLogPts, alpha, c, xFitVals, yFitVals, fitLine,
   xRange, yRange, epilog,
   idxFMin, idxFMax, xLineMin, xLineMax, yLogMin, yLogMax, vertLinePrims},

  th = ThemeFromOpts[{opts}];
  fg = th["Fg"];

  freq  = stats["FreqHz"];
  yMean = stats["YMean"];
  yVar  = stats["YVar"];

  logFreq = Log10[freq];
  logFMin = Min[logFreq]; logFMax = Max[logFreq];

  goodMask  = MapThread[#1 > 0 && #2 > 0 &, {yMean, yVar}];
  yMeanGood = Pick[yMean, goodMask];
  yVarGood  = Pick[yVar,  goodMask];
  freqGood  = Pick[freq,  goodMask];

  colsGood        = ColorData["Rainbow"] /@ Rescale[Log10[freqGood], {logFMin, logFMax}];
  singletonSets   = {{#1, #2}} & @@@ Transpose[{yMeanGood, yVarGood}];
  singletonStyles = Directive[#, PointSize[0.013]] & /@ colsGood;

  fitMask  = (# >= fMin && # <= fMax) & /@ freqGood;
  yMeanFit = Pick[yMeanGood, fitMask];
  yVarFit  = Pick[yVarGood,  fitMask];

  goodPairs  = Select[Transpose[{yMeanFit, yVarFit}], #[[1]] > 0 && #[[2]] > 0 &];
  goodLogPts = {Log10[#[[1]]], Log10[#[[2]]]} & /@ goodPairs;

  {alpha, c} = If[Length[goodLogPts] >= 3,
    Module[{lmFit, p},
      lmFit = LinearModelFit[goodLogPts, {1, x}, x];
      p = lmFit["BestFitParameters"];
      {p[[2]], 10.^p[[1]]}
    ],
    {2., 1.}
  ];

  xFitVals = 10.^Range[Log10[Min[yMeanFit]] - 0.05, Log10[Max[yMeanFit]] + 0.05, 0.02];
  yFitVals = c * xFitVals^alpha;
  fitLine  = Transpose[{xFitVals, yFitVals}];

  xRange = MinMax[yMeanGood];
  yRange = MinMax[yVarGood];

  (* vertical dashed lines at fMin/fMax — Log10 (scaled-axis) coordinates *)
  idxFMin  = First @ Ordering[Abs[Log10[freqGood] - Log10[fMin]], 1];
  idxFMax  = First @ Ordering[Abs[Log10[freqGood] - Log10[fMax]], 1];
  xLineMin = Log10[yMeanGood[[idxFMin]]];
  xLineMax = Log10[yMeanGood[[idxFMax]]];
  yLogMin  = Log10[Min[yVarGood]];
  yLogMax  = Log10[Max[yVarGood]];
  vertLinePrims = {
    Dashed, GrayLevel[0.65], AbsoluteThickness[1.5],
    Line[{{xLineMin, yLogMin}, {xLineMin, yLogMax}}],
    Line[{{xLineMax, yLogMin}, {xLineMax, yLogMax}}]
  };

  epilog = Join[
    {Style[Text[annotLabelFor[alpha], Scaled[{0.04, 0.96}], {-1, 1}], fg, 12]},
    {vertLinePrims},
    freqLabelsFor[freqGood, yMeanGood, yVarGood, fg]
  ];

  variancePowerLawShow[th, singletonSets, singletonStyles, fitLine, xRange, yRange, epilog]
]

End[]
EndPackage[]
