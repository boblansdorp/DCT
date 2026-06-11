(* ::Package:: *)

(* DCTPlots`PlotTuning
   Parameter-tuning diagnostic plots (Section 3 of the analysis notebook).

   Public API:
     PlotLCurve[lCurveData]            L-curve: roughness (log x) vs weighted
     PlotLCurve[lCurveData, labels]    misfit (linear y); optional per-point labels.

   Theme-aware: default DCTPlots`$DCTTheme ("Dark" | "Publication"); override
   one call with Theme -> "Publication".
*)

BeginPackage["DCTPlots`"]

PlotLCurve::usage =
  "PlotLCurve[lCurveData] plots the regularization L-curve: roughness \
||D^2 g|| (log x) vs weighted admittance misfit (linear y), one point per \
lambda. PlotLCurve[lCurveData, labels] annotates each point with a label \
(e.g. file name); labels = None suppresses them.";

PlotLambdaResiduals::usage =
  "PlotLambdaResiduals[lSweepSpecs] returns a stacked pair of Bode-residual \
panels (|dY|/|Y| in %, and dPhase in degrees, vs frequency), one rainbow trace \
per lambda in the sweep \[Dash] for verifying the chosen regularization lambda.";

Begin["`Private`"]

Options[PlotLCurve] = {Theme -> Automatic};

PlotLCurve[lCurveData_List, labels : (_List | None) : None,
           opts : OptionsPattern[]] := Module[
  {th, fg, plotOpts, xR, yR, xLo, xHi, yt0, ymaj, dy, yLo, yHi,
   pointLabels, epilog},

  th       = ThemeFromOpts[{opts}];
  fg       = th["Fg"];
  plotOpts = FilterRules[{opts}, Options[ListLinePlot]];

  xR = MinMax[lCurveData[[All, 1]]];   (* roughness — log x  *)
  yR = MinMax[lCurveData[[All, 2]]];   (* misfit    — linear y *)

  (* Frame bounds = the ticks' own extent so ticks/numbers reach the frame
     edges: x to the enclosing decades; y out by one major-tick step past
     the data (the data then sits comfortably inside). *)
  xLo  = 10.^Floor[Log10[xR[[1]]]];
  xHi  = 10.^Ceiling[Log10[xR[[2]]]];
  yt0  = ThemeLinTicks[yR[[1]], yR[[2]], fg];
  ymaj = Sort @ Cases[yt0, {p_, l_, ___} /; l =!= "" :> p];
  dy   = If[Length[ymaj] >= 2, ymaj[[-1]] - ymaj[[-2]], yR[[2]] - yR[[1]]];
  yLo  = If[ymaj === {}, yR[[1]], ymaj[[1]]  - dy];
  yHi  = If[ymaj === {}, yR[[2]], ymaj[[-1]] + dy];

  (* point labels placed at Scaled fractions of the (log-x, linear-y) frame.
     Scaled coords render reliably on the native log plot, where data-coordinate
     Epilog primitives silently vanish at this tiny y-scale (~1e-4). *)
  pointLabels = If[labels === None, {},
    MapThread[
      Style[Text[FileNameTake[ToString[#2]],
        Scaled[{(Log10[#1[[1]]] - Log10[xLo]) / (Log10[xHi] - Log10[xLo]),
                (#1[[2]] - yLo) / (yHi - yLo)}], {-1.3, 0.5}], fg, 10] &,
      {lCurveData, labels}]
  ];

  epilog = Join[
    pointLabels,
    {Style[Text["\[LeftArrow] smoother",   Scaled[{0.12, 0.95}], {-1, 1}], Gray, 11],
     Style[Text["better fit \[RightArrow]", Scaled[{0.88, 0.05}], { 1, -1}], Gray, 11]}
  ];

  ListLogLinearPlot[
    lCurveData,
    plotOpts,
    Joined           -> True,
    PlotStyle        -> Directive[fg, Opacity[0.6], AbsoluteThickness[1.5]],
    PlotMarkers      -> {Graphics[{fg, Disk[{0, 0}, 1]}], 0.03},
    Sequence @@ ThemeChrome[th, 13, 1.2],
    FrameTicks       -> {{ThemeLinTicks[yLo, yHi, fg], None},
                         {ThemeLogTicks[xLo, xHi, fg], None}},
    FrameTicksStyle  -> Directive[fg, 12],
    FrameLabel       -> {
      Style["Roughness  \[LeftDoubleBracketingBar]D\[CenterDot]D\[CenterDot]g\[RightDoubleBracketingBar]  (F/dec)", 14, fg],
      Style["Weighted misfit  \[LeftDoubleBracketingBar]\[Sqrt]w \[CenterDot] \[CapitalDelta]Y\[RightDoubleBracketingBar]", 14, fg]
    },
    Epilog           -> epilog,
    PlotRange        -> {{xLo, xHi}, {yLo, yHi}},
    PlotRangePadding -> None,
    ImageSize        -> 500,
    AspectRatio      -> 0.8
  ]
]

(* ------------------------------------------------------------------ *)

Options[PlotLambdaResiduals] = {Theme -> Automatic};

PlotLambdaResiduals[lSweepSpecs_List, opts : OptionsPattern[]] := Module[
  {th, fg, plotOpts, n, cols, styles, magTraces, phaseTraces, freq,
   fMin, fMax, magMax, phMin, phMax, panel},

  th       = ThemeFromOpts[{opts}];
  fg       = th["Fg"];
  plotOpts = FilterRules[{opts}, Options[ListLinePlot]];

  n      = Length[lSweepSpecs];
  cols   = Table[ColorData["Rainbow"][u], {u, 0.05, 0.95, 0.90 / Max[n - 1, 1]}];
  styles = Directive[#, AbsoluteThickness[2.]] & /@ cols;

  magTraces = Table[
    Module[{sp = lSweepSpecs[[i, "Spec"]]},
      Transpose[{sp["FreqHz"],
        100. Abs[sp["YIntFit"] - sp["YIntData"]] / Abs[sp["YIntData"]]}]],
    {i, n}];
  phaseTraces = Table[
    Module[{sp = lSweepSpecs[[i, "Spec"]]},
      Transpose[{sp["FreqHz"],
        (Arg[sp["YIntFit"]] - Arg[sp["YIntData"]]) * 180/Pi}]],
    {i, n}];

  freq         = lSweepSpecs[[1, "Spec", "FreqHz"]];
  {fMin, fMax} = MinMax[freq];
  magMax       = Max[Flatten[magTraces[[All, All, 2]]]] * 1.05;
  {phMin, phMax} = MinMax[Flatten[phaseTraces[[All, All, 2]]]];

  panel[traces_, yLabel_, yLo_, yHi_] := ListLogLinearPlot[
    traces,
    plotOpts,
    Joined           -> True,
    PlotStyle        -> styles,
    Sequence @@ ThemeChrome[th, 13, 1.2],
    FrameTicks       -> {{ThemeLinTicks[yLo, yHi, fg], None},
                         {ThemeLogTicks[fMin, fMax, fg], None}},
    FrameTicksStyle  -> Directive[fg, 12],
    FrameLabel       -> {Style["f (Hz)", 13, fg], Style[yLabel, 13, fg]},
    PlotRange        -> {Automatic, {yLo, yHi}},
    ImageSize        -> 360,
    AspectRatio      -> 0.6];

  Column[{
    panel[magTraces,   "\[LeftBracketingBar]\[CapitalDelta]Y\[RightBracketingBar] / \[LeftBracketingBar]Y\[RightBracketingBar] (%)", 0., magMax],
    panel[phaseTraces, "\[CapitalDelta]Phase (\[Degree])", phMin, phMax]
  }, Spacings -> 0.3]
]

End[]
EndPackage[]
