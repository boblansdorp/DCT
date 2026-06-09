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

(* Shared style constants (same as rest of DCTPlots) *)
$darkBg2   = GrayLevel[0.12];
$fgStyle2  = Directive[White, AbsoluteThickness[1.2]];
$lblStyle2 = Directive[White, 14, FontFamily -> "Arial"];

logTicksFromList[vals_List] := Module[{lo, hi},
  lo = Floor[Log10[Min[Select[vals, Positive]]]];
  hi = Ceiling[Log10[Max[Select[vals, Positive]]]];
  Table[{10.^p, Superscript[10, p]}, {p, lo, hi}]
]

(* ------------------------------------------------------------------ *)

PlotCV[stats_Association, cvThresh_: 0.1] := Module[
  {freq, cv, cvPairs, threshPair, xTicks, yMax},

  freq = stats["FreqHz"];
  cv   = stats["CV"];

  cvPairs    = Transpose[{freq, cv}];
  threshPair = {{Min[freq], cvThresh}, {Max[freq], cvThresh}};
  xTicks     = logTicksFromList[freq];
  yMax       = Max[Select[cv, NumericQ]] * 1.1;

  ListLinePlot[
    {cvPairs, threshPair},
    ScalingFunctions  -> {"Log10", None},
    PlotStyle         -> {
      Directive[ColorData["Rainbow"][0.55], AbsoluteThickness[2.2]],
      Directive[Red, Dashed, AbsoluteThickness[1.8]]
    },
    Frame             -> True,
    Axes              -> False,
    Background        -> $darkBg2,
    FrameStyle        -> $fgStyle2,
    FrameTicks        -> {{Automatic, None}, {xTicks, None}},
    FrameTicksStyle   -> Directive[White, 13],
    LabelStyle        -> $lblStyle2,
    FrameLabel        -> {
      Style["Frequency (Hz)", 14, White],
      Style["\[Sigma]/\[Mu] of |Y(f)|", 14, White]
    },
    PlotRange         -> {{Min[freq], Max[freq]}, {0, yMax}},
    PlotRangePadding  -> {{Scaled[0.02], Scaled[0.02]}, {0, Scaled[0.05]}},
    ImageSize         -> 600,
    AspectRatio       -> 0.55,
    PlotLegends       -> Placed[
      LineLegend[
        {Directive[ColorData["Rainbow"][0.55], AbsoluteThickness[2.2]],
         Directive[Red, Dashed, AbsoluteThickness[1.8]]},
        {"CV(f)", "Threshold " <> ToString[cvThresh]},
        LabelStyle      -> Directive[White, 12],
        Background      -> $darkBg2,
        LegendFunction  -> (Framed[#, Background -> $darkBg2,
                              FrameStyle -> $fgStyle2] &)
      ],
      {Right, Top}
    ]
  ]
]

(* ------------------------------------------------------------------ *)

PlotVariancePowerLaw[stats_Association] := Module[
  {freq, yMean, yVar, alpha, c, logFreq, logFMin, logFMax,
   goodMask, yMeanGood, yVarGood, freqGood,
   colsGood, singletonSets, singletonStyles,
   xFitVals, yFitVals, fitLine,
   xTicks, yTicks, annotLabel,
   freqTargets, freqLabelPts, freqLabelPrims},

  freq  = stats["FreqHz"];
  yMean = stats["YMean"];
  yVar  = stats["YVar"];
  alpha = stats["PowerLawAlpha"];
  c     = stats["PowerLawC"];

  logFreq = Log10[freq];
  logFMin = Min[logFreq]; logFMax = Max[logFreq];

  (* keep only finite positive points *)
  goodMask  = MapThread[#1 > 0 && #2 > 0 &, {yMean, yVar}];
  yMeanGood = Pick[yMean, goodMask];
  yVarGood  = Pick[yVar,  goodMask];
  freqGood  = Pick[freq,  goodMask];

  (* per-point colour by log10(f) *)
  colsGood = ColorData["Rainbow"] /@
    Rescale[Log10[freqGood], {logFMin, logFMax}];

  (* N single-point datasets for per-point colouring *)
  singletonSets   = {{#1, #2}} & @@@ Transpose[{yMeanGood, yVarGood}];
  singletonStyles = Directive[#, PointSize[0.013]] & /@ colsGood;

  (* fit line spanning data range with a small margin *)
  xFitVals = 10.^Range[
    Log10[Min[yMeanGood]] - 0.05,
    Log10[Max[yMeanGood]] + 0.05,
    0.02
  ];
  yFitVals = c * xFitVals^alpha;
  fitLine  = Transpose[{xFitVals, yFitVals}];

  xTicks = logTicksFromList[yMeanGood];
  yTicks = logTicksFromList[yVarGood];

  (* frequency labels at decade marks — only include if a data point is within 0.3 decades *)
  freqTargets = {1., 10., 100., 1000., 10000.};
  freqLabelPts = Map[
    Function[f,
      Module[{idx, logDiff},
        idx     = First @ Ordering[Abs[Log10[freqGood] - Log10[f]], 1];
        logDiff = Abs[Log10[freqGood[[idx]]] - Log10[f]];
        If[logDiff < 0.3,
          {yMeanGood[[idx]], yVarGood[[idx]], f},
          Nothing
        ]
      ]
    ],
    freqTargets
  ];
  freqLabelPrims = Style[
    Text[Row[{ToString[Round[#[[3]]]], " Hz"}], {Log[#[[1]]], Log[#[[2]]]}, {-1.3, 0.6}],
    White, 10
  ] & /@ freqLabelPts;

  annotLabel = Row[{
    "Var \[Proportional] |Y|",
    Superscript["\[Alpha]", ""],
    ",  \[Alpha] = ", NumberForm[alpha, {4, 2}],
    "   \[DoubleRightArrow]   weightPower = ",
    NumberForm[-alpha, {4, 2}]
  }];

  Show[
    ListLogLogPlot[
      singletonSets,
      PlotStyle         -> singletonStyles,
      Frame             -> True,
      Axes              -> False,
      Background        -> $darkBg2,
      FrameStyle        -> $fgStyle2,
      FrameTicks        -> {{yTicks, None}, {xTicks, None}},
      FrameTicksStyle   -> Directive[White, 13],
      LabelStyle        -> $lblStyle2,
      FrameLabel        -> {
        {Style["Var(|Y(f)|) (S^2)", 14, White], None},
        {Style["\[LeftAngleBracket]|Y(f)|\[RightAngleBracket] (S)", 14, White], None}
      },
      PlotRangePadding  -> Scaled[0.04],
      ImageSize         -> 600,
      AspectRatio       -> 0.65
    ],
    ListLogLogPlot[
      {fitLine},
      PlotStyle -> Directive[White, AbsoluteThickness[2.2]],
      Joined    -> True
    ],
    Epilog -> Join[
      {Style[Text[annotLabel, Scaled[{0.04, 0.96}], {-1, 1}], White, 12]},
      freqLabelPrims
    ]
  ]
]

(* ------------------------------------------------------------------ *)

PlotVariancePowerLaw[stats_Association, fMin_?NumericQ, fMax_?NumericQ] := Module[
  {freq, yMean, yVar, logFreq, logFMin, logFMax,
   goodMask, yMeanGood, yVarGood, freqGood,
   fitMask, yMeanFit, yVarFit,
   colsGood, singletonSets, singletonStyles,
   goodPairs, goodLogPts, alpha, c,
   xFitVals, yFitVals, fitLine,
   xTicks, yTicks, annotLabel,
   freqTargets, freqLabelPts, freqLabelPrims,
   idxFMin, idxFMax, xLineMin, xLineMax, yLogMin, yLogMax, vertLinePrims},

  freq  = stats["FreqHz"];
  yMean = stats["YMean"];
  yVar  = stats["YVar"];

  logFreq = Log10[freq];
  logFMin = Min[logFreq]; logFMax = Max[logFreq];

  (* all clean points for scatter *)
  goodMask  = MapThread[#1 > 0 && #2 > 0 &, {yMean, yVar}];
  yMeanGood = Pick[yMean, goodMask];
  yVarGood  = Pick[yVar,  goodMask];
  freqGood  = Pick[freq,  goodMask];

  colsGood        = ColorData["Rainbow"] /@ Rescale[Log10[freqGood], {logFMin, logFMax}];
  singletonSets   = {{#1, #2}} & @@@ Transpose[{yMeanGood, yVarGood}];
  singletonStyles = Directive[#, PointSize[0.013]] & /@ colsGood;

  (* power-law fit on [fMin, fMax] range only *)
  fitMask   = (# >= fMin && # <= fMax) & /@ freqGood;
  yMeanFit  = Pick[yMeanGood, fitMask];
  yVarFit   = Pick[yVarGood,  fitMask];

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

  xFitVals = 10.^Range[
    Log10[Min[yMeanFit]] - 0.05,
    Log10[Max[yMeanFit]] + 0.05,
    0.02
  ];
  yFitVals = c * xFitVals^alpha;
  fitLine  = Transpose[{xFitVals, yFitVals}];

  xTicks = logTicksFromList[yMeanGood];
  yTicks = logTicksFromList[yVarGood];

  (* vertical dashed lines at fMin/fMax — in Log[] (natural log) coordinates *)
  idxFMin  = First @ Ordering[Abs[Log10[freqGood] - Log10[fMin]], 1];
  idxFMax  = First @ Ordering[Abs[Log10[freqGood] - Log10[fMax]], 1];
  xLineMin = Log[yMeanGood[[idxFMin]]];
  xLineMax = Log[yMeanGood[[idxFMax]]];
  yLogMin  = Log[Min[yVarGood]];
  yLogMax  = Log[Max[yVarGood]];
  vertLinePrims = {
    Dashed, GrayLevel[0.65], AbsoluteThickness[1.5],
    Line[{{xLineMin, yLogMin}, {xLineMin, yLogMax}}],
    Line[{{xLineMax, yLogMin}, {xLineMax, yLogMax}}]
  };

  freqTargets = {1., 10., 100., 1000., 10000.};
  freqLabelPts = Map[
    Function[f,
      Module[{idx, logDiff},
        idx     = First @ Ordering[Abs[Log10[freqGood] - Log10[f]], 1];
        logDiff = Abs[Log10[freqGood[[idx]]] - Log10[f]];
        If[logDiff < 0.3, {yMeanGood[[idx]], yVarGood[[idx]], f}, Nothing]
      ]
    ],
    freqTargets
  ];
  freqLabelPrims = Style[
    Text[Row[{ToString[Round[#[[3]]]], " Hz"}], {Log[#[[1]]], Log[#[[2]]]}, {-1.3, 0.6}],
    White, 10
  ] & /@ freqLabelPts;

  annotLabel = Row[{
    "Var \[Proportional] |Y|",
    Superscript["\[Alpha]", ""],
    ",  \[Alpha] = ", NumberForm[alpha, {4, 2}],
    "   \[DoubleRightArrow]   weightPower = ",
    NumberForm[-alpha, {4, 2}]
  }];

  Show[
    ListLogLogPlot[
      singletonSets,
      PlotStyle         -> singletonStyles,
      Frame             -> True,
      Axes              -> False,
      Background        -> $darkBg2,
      FrameStyle        -> $fgStyle2,
      FrameTicks        -> {{yTicks, None}, {xTicks, None}},
      FrameTicksStyle   -> Directive[White, 13],
      LabelStyle        -> $lblStyle2,
      FrameLabel        -> {
        {Style["Var(|Y(f)|) (S^2)", 14, White], None},
        {Style["\[LeftAngleBracket]|Y(f)|\[RightAngleBracket] (S)", 14, White], None}
      },
      PlotRangePadding  -> Scaled[0.04],
      ImageSize         -> 600,
      AspectRatio       -> 0.65
    ],
    ListLogLogPlot[
      {fitLine},
      PlotStyle -> Directive[White, AbsoluteThickness[2.2]],
      Joined    -> True
    ],
    Epilog -> Join[
      {Style[Text[annotLabel, Scaled[{0.04, 0.96}], {-1, 1}], White, 12]},
      {vertLinePrims},
      freqLabelPrims
    ]
  ]
]

End[]
EndPackage[]
