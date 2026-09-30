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

  xR = MinMax[lCurveData[[All, 1]]];   (* roughness, log x  *)
  yR = MinMax[lCurveData[[All, 2]]];   (* misfit,    linear y *)

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
  (* labels may be strings (file names, shortened) or typeset expressions such as
     Row[{Subscript["\[Lambda]", "ND"], " = 0.3"}], passed through. Labels sit to the
     right of their point, except in the flat tail (right third of the frame) where
     the points crowd together: those are stacked in rows above-left of their points,
     one row higher for each successive tail point, in the free space above the curve. *)
  pointLabels = If[labels === None, {},
    Module[{sxs, sys, lbl, tail, yBase, tailIdx = 0},
      sxs = (Log10[lCurveData[[All, 1]]] - Log10[xLo]) / (Log10[xHi] - Log10[xLo]);
      sys = (lCurveData[[All, 2]] - yLo) / (yHi - yLo);
      lbl = If[StringQ[#], FileNameTake[#], #] & /@ labels;
      tail  = Select[Range[Length[sxs]], sxs[[#]] > 0.62 &];
      (* tail labels sit in rows starting a little above the highest tail point *)
      yBase = If[tail === {}, 0., 0.12 + Max[sys[[tail]]]];
      Table[
        If[sxs[[i]] <= 0.62,
          (* steep part: label to the right of the point. Around the knee (right half,
             lower third) the right side is where the tail leaders rise, so those two
             or three labels go to the left, into the empty concave side. *)
          If[sxs[[i]] > 0.5 && sys[[i]] < 0.35,
            Style[Text[lbl[[i]], Scaled[{sxs[[i]], sys[[i]]}], {1.3, 0.5}], fg, 10],
            (* to the right and a little above, clear of the marker and of the curve
               continuing down and to the right *)
            Style[Text[lbl[[i]], Scaled[{sxs[[i]], sys[[i]]}], {-1.5, -0.7}], fg, 10]],
          With[{y = yBase + 0.07 (tailIdx++)},
            {{GrayLevel[0.6], AbsoluteThickness[0.6],
              Line[{Scaled[{sxs[[i]], sys[[i]] + 0.03}], Scaled[{sxs[[i]], y - 0.025}]}]},
             (* centred above the point, or right-aligned at the frame edge if too far right *)
             If[sxs[[i]] > 0.88,
               Style[Text[lbl[[i]], Scaled[{0.99, y}], {1, 0}], fg, 10],
               Style[Text[lbl[[i]], Scaled[{sxs[[i]], y}], {0, 0}], fg, 10]]}]],
        {i, Length[lCurveData]}]]
  ];

  epilog = Join[
    pointLabels,
    (* left margin below the top point: clear of that point's label *)
    {Style[Text["\[LeftArrow] smoother",   Scaled[{0.02, 0.35}], {-1, 0}], Gray, 11],
     (* misfit is the y axis, so better fit is DOWN, not right *)
     Style[Text["better fit \[DownArrow]", Scaled[{0.88, 0.05}], { 1, -1}], Gray, 11]}
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
      (* D denotes the second-difference matrix (rows 1, -2, 1), applied once, matching
         the manuscript notation and the Figure S4 caption: roughness = ||D.g||. The dot
         is the matrix-vector product; the bare stencil carries no 1/Delta^2, so the
         units are those of g, F/decade. *)
      Style["Roughness  \[LeftDoubleBracketingBar]D\[CenterDot]g\[RightDoubleBracketingBar]  (F/decade)", 14, fg],
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
