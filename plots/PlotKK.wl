(* ::Package:: *)

(* DCTPlots`PlotKK
   Visualisation of linearized Kramers-Kronig residuals from KKCheck[].

   Public API:
     PlotKKResiduals[kkResult]             residuals for one spectrum
     PlotKKResiduals[kkResult, threshold]  same, custom threshold bands

   Theme-aware: default DCTPlots`$DCTTheme ("Dark" | "Publication"); override
   one call with Theme -> "Publication". The Re/Im traces keep their fixed
   blue/orange colours; chrome, threshold band and legend follow the theme.
*)

BeginPackage["DCTPlots`"]

PlotKKResiduals::usage =
  "PlotKKResiduals[kkResult, threshold:0.02] plots the relative admittance \
residuals (real and imaginary) vs frequency from KKCheck[]. Dashed lines mark \
the pass/fail threshold."

PlotKKFMinHistogram::usage =
  "PlotKKFMinHistogram[kkMins] returns a histogram of the per-file linearized \
Kramers-Kronig valid lower-bound frequencies (a long right tail flags outlier files)."

Begin["`Private`"]

$kkReStyle = Directive[RGBColor[0.35, 0.75, 1.],  AbsoluteThickness[1.8]];  (* data *)
$kkImStyle = Directive[RGBColor[1.,   0.65, 0.15], AbsoluteThickness[1.8]];  (* data *)

Options[PlotKKResiduals] = {Theme -> Automatic};

PlotKKResiduals[kkResult_Association, threshold : (_?NumericQ) : 0.02,
                opts : OptionsPattern[]] := Module[
  {th, fg, plotOpts, freq, resRe, resIm, ptsRe, ptsIm, fMin, fMax, yMax, bandStyle},

  th       = ThemeFromOpts[{opts}];
  fg       = th["Fg"];
  plotOpts = FilterRules[{opts}, Options[ListLinePlot]];

  freq   = kkResult["FreqHz"];
  resRe  = kkResult["YResRe"];
  resIm  = kkResult["YResIm"];
  ptsRe  = Transpose[{freq, resRe}];
  ptsIm  = Transpose[{freq, resIm}];
  {fMin, fMax} = MinMax[freq];
  yMax   = Max[1.5 threshold, 1.2 Max[Abs[resRe], Abs[resIm]]];

  (* threshold band drawn as 2-point data series so ScalingFunctions positions
     it (real frequency coords) — avoids the old Show + Log[] coordinate hack
     and lets the theme outward-tick generators apply. *)
  bandStyle = Directive[fg, Opacity[0.35], Dashing[{0.015, 0.01}], AbsoluteThickness[1.]];

  ListLinePlot[
    {ptsRe, ptsIm,
     {{fMin, threshold},  {fMax, threshold}},
     {{fMin, -threshold}, {fMax, -threshold}}},
    plotOpts,
    ScalingFunctions -> {"Log10", None},
    PlotStyle        -> {$kkReStyle, $kkImStyle, bandStyle, bandStyle},
    Sequence @@ ThemeChrome[th, 13, 1.1],
    FrameTicks       -> {{ThemeLinTicks[-yMax, yMax, fg], None},
                         {ThemeLogTicks[fMin, fMax, fg], None}},
    FrameTicksStyle  -> Directive[fg, 12],
    FrameLabel       -> {
      Style["Frequency (Hz)", 13, fg],
      Style["Relative admittance residual", 13, fg]
    },
    PlotRange   -> {All, {-yMax, yMax}},
    GridLines   -> {None, {0}},
    GridLinesStyle -> Directive[fg, Opacity[0.2], AbsoluteThickness[0.7]],
    AspectRatio -> 0.65,
    ImageSize   -> 620,
    PlotLegends -> Placed[LineLegend[
      {$kkReStyle, $kkImStyle},
      {"\[Delta]Re(Y)/|Y|", "\[Delta]Im(Y)/|Y|"},
      LegendMarkerSize -> 20,
      LabelStyle       -> Directive[fg, 12],
      Background       -> th["Bg"],
      LegendFunction   -> (Framed[#, Background -> th["Bg"],
        FrameStyle -> Directive[fg, AbsoluteThickness[1.1]]] &)],
      {0.82, 0.88}]
  ]
]

(* ------------------------------------------------------------------ *)

Options[PlotKKFMinHistogram] = {Theme -> Automatic};

PlotKKFMinHistogram[kkMins_List, opts : OptionsPattern[]] := Module[
  {th, fg, plotOpts, edges, counts, xMin, xMax, yMax},
  th       = ThemeFromOpts[{opts}];
  fg       = th["Fg"];
  plotOpts = FilterRules[{opts}, Options[Histogram]];

  (* bin counts/edges (same Automatic spec as the Histogram) give the tick
     ranges so the outward generators can be applied to a Histogram *)
  {edges, counts} = HistogramList[kkMins, Automatic, "Count"];
  xMin = First[edges];  xMax = Last[edges];
  yMax = Max[counts];

  Histogram[kkMins, Automatic, "Count",
    plotOpts,
    Sequence @@ ThemeChrome[th, 13, 1.1],
    ChartStyle      -> Directive[ColorData["Rainbow"][0.55], EdgeForm[Directive[fg, AbsoluteThickness[0.6]]]],
    FrameTicks      -> {{ThemeLinTicks[0., yMax, fg], None},
                        {ThemeLinTicks[xMin, xMax, fg], None}},
    FrameTicksStyle -> Directive[fg, 12],
    FrameLabel      -> {Style["Lin-KK fMin (Hz)", 13, fg], Style["Files", 13, fg]},
    PlotLabel       -> Style["Valid lower bound \[LongDash] all drift files", fg, 13],
    ImageSize       -> 440
  ]
]

End[]
EndPackage[]
