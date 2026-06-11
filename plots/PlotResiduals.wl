(* ::Package:: *)

(* DCTPlots`PlotResiduals
   Residual analysis plots comparing raw data to the DCT fit.

   Public API:
     PlotResidualMag[spec]     |(Z_data - Z_fit) / Z_data| vs frequency
     PlotResidualPhase[spec]   |phase_data - phase_fit| in degrees vs f
     PlotResiduals[spec]       GraphicsGrid: magnitude + phase side by side

   Theme-aware: default DCTPlots`$DCTTheme ("Dark" | "Publication"); override
   one call with Theme -> "Publication". See plots/PlotTheme.wl.
*)

BeginPackage["DCTPlots`"]

PlotResidualMag::usage =
  "PlotResidualMag[spec] plots the relative magnitude residual \
|(Z_data - Z_fit)/Z_data| vs frequency on a log-log scale."

PlotResidualPhase::usage =
  "PlotResidualPhase[spec] plots the absolute phase residual \
|phase_data - phase_fit| (degrees) vs frequency on a semilog scale."

PlotResiduals::usage =
  "PlotResiduals[spec] returns a GraphicsGrid with magnitude and phase \
residual panels side by side."

Begin["`Private`"]

resStyle     = Directive[Lighter[Blue, 0.4], AbsoluteThickness[1.8]];  (* data colour *)
resImageSize = 420;
resAspect    = 0.72;

Options[PlotResidualMag] = {Theme -> Automatic};

PlotResidualMag[spec_Association, opts : OptionsPattern[]] := Module[
  {th, fg, plotOpts, freq, zData, zFit, pts, fMin, fMax},

  th       = ThemeFromOpts[{opts}];
  fg       = th["Fg"];
  plotOpts = FilterRules[{opts}, Options[ListLinePlot]];

  freq  = spec["FreqHz"];
  zData = spec["ZData"];
  zFit  = spec["ZFit"];
  pts   = Transpose[{freq, Abs[(zData - zFit) / zData]}];
  {fMin, fMax} = MinMax[freq];

  ListLogLogPlot[
    pts,
    plotOpts,
    Joined           -> True,
    PlotStyle        -> resStyle,
    Sequence @@ ThemeChrome[th, 14, 1.1],
    FrameTicks       -> {{ThemeLogTicks[10.^-4, 1., fg], None},
                         {ThemeLogTicks[fMin, fMax, fg], None}},
    FrameTicksStyle  -> Directive[fg, 12],
    FrameLabel       -> {
      Style["Frequency (Hz)", 14, fg],
      Style["|(Z_data - Z_fit)/Z_data|", 14, fg]
    },
    PlotRange   -> {All, {10^-4, 1}},
    AspectRatio -> resAspect,
    ImageSize   -> resImageSize
  ]
]

Options[PlotResidualPhase] = {Theme -> Automatic};

PlotResidualPhase[spec_Association, opts : OptionsPattern[]] := Module[
  {th, fg, plotOpts, freq, zData, zFit, pts, fMin, fMax, yMax},

  th       = ThemeFromOpts[{opts}];
  fg       = th["Fg"];
  plotOpts = FilterRules[{opts}, Options[ListLinePlot]];

  freq  = spec["FreqHz"];
  zData = spec["ZData"];
  zFit  = spec["ZFit"];
  pts   = Transpose[{
    freq,
    Abs[(180. / Pi) Arg /@ zData - (180. / Pi) Arg /@ zFit]
  }];
  {fMin, fMax} = MinMax[freq];
  yMax = Max[pts[[All, 2]]] * 1.05;

  ListLogLinearPlot[
    pts,
    plotOpts,
    Joined           -> True,
    PlotStyle        -> resStyle,
    Sequence @@ ThemeChrome[th, 14, 1.1],
    FrameTicks       -> {{ThemeLinTicks[0., yMax, fg], None},
                         {ThemeLogTicks[fMin, fMax, fg], None}},
    FrameTicksStyle  -> Directive[fg, 12],
    FrameLabel       -> {
      Style["Frequency (Hz)", 14, fg],
      Style["|phase residual| (\[Degree])", 14, fg]
    },
    PlotRange   -> {Automatic, {0, yMax}},
    AspectRatio -> resAspect,
    ImageSize   -> resImageSize
  ]
]

Options[PlotResiduals] = {Theme -> Automatic};

PlotResiduals[spec_Association, opts : OptionsPattern[]] := Module[
  {th},
  th = ThemeFromOpts[{opts}];
  GraphicsGrid[
    {{PlotResidualMag[spec, Theme -> th], PlotResidualPhase[spec, Theme -> th]}},
    ImageSize  -> 900,
    Spacings   -> {0.5, 0.5},
    Background -> th["Bg"]
  ]
]

(* ------------------------------------------------------------------ *)
(* Multi-spectrum overlays: every fit at once, semi-transparent, with  *)
(* reference lines at 1% magnitude / 1 degree phase (the paper claim). *)
(* ------------------------------------------------------------------ *)

resOverlayStyle = Directive[Lighter[Blue, 0.2], Opacity[0.2], AbsoluteThickness[0.9]];
resRefStyle[fg_] := Directive[fg, Dashing[{0.018, 0.014}], AbsoluteThickness[1.4]];

PlotResidualMag[specs : {__Association}, opts : OptionsPattern[]] := Module[
  {th, fg, plotOpts, traces, fMin, fMax},

  th       = ThemeFromOpts[{opts}];
  fg       = th["Fg"];
  plotOpts = FilterRules[{opts}, Options[ListLinePlot]];

  traces = Function[s,
    Transpose[{s["FreqHz"], Abs[(s["ZData"] - s["ZFit"]) / s["ZData"]]}]] /@ specs;
  {fMin, fMax} = MinMax[Flatten[traces[[All, All, 1]]]];

  ListLogLogPlot[
    traces,
    plotOpts,
    Joined           -> True,
    PlotStyle        -> resOverlayStyle,
    Sequence @@ ThemeChrome[th, 14, 1.1],
    GridLines        -> {None, {{0.05, resRefStyle[fg]}}},   (* 5% amplitude bound *)
    FrameTicks       -> {{ThemeLogTicks[10.^-4, 1., fg], None},
                         {ThemeLogTicks[fMin, fMax, fg], None}},
    FrameTicksStyle  -> Directive[fg, 12],
    FrameLabel       -> {
      Style["Frequency (Hz)", 14, fg],
      Style["|(Z_data - Z_fit)/Z_data|", 14, fg]
    },
    PlotRange   -> {All, {10^-4, 1}},
    AspectRatio -> resAspect,
    ImageSize   -> resImageSize
  ]
]

PlotResidualPhase[specs : {__Association}, opts : OptionsPattern[]] := Module[
  {th, fg, plotOpts, traces, fMin, fMax, yMax},

  th       = ThemeFromOpts[{opts}];
  fg       = th["Fg"];
  plotOpts = FilterRules[{opts}, Options[ListLinePlot]];

  traces = Function[s,
    Transpose[{s["FreqHz"],
      Abs[(180. / Pi) Arg /@ s["ZData"] - (180. / Pi) Arg /@ s["ZFit"]]}]] /@ specs;
  {fMin, fMax} = MinMax[Flatten[traces[[All, All, 1]]]];
  (* headroom so the 3 degree reference line is always visible *)
  yMax = Min[6., Max[3.4, Max[Flatten[traces[[All, All, 2]]]] * 1.1]];

  ListLogLinearPlot[
    traces,
    plotOpts,
    Joined           -> True,
    PlotStyle        -> resOverlayStyle,
    Sequence @@ ThemeChrome[th, 14, 1.1],
    GridLines        -> {None, {{3., resRefStyle[fg]}}},   (* 3 degree phase bound *)
    FrameTicks       -> {{ThemeLinTicks[0., yMax, fg], None},
                         {ThemeLogTicks[fMin, fMax, fg], None}},
    FrameTicksStyle  -> Directive[fg, 12],
    FrameLabel       -> {
      Style["Frequency (Hz)", 14, fg],
      Style["|phase residual| (\[Degree])", 14, fg]
    },
    PlotRange   -> {Automatic, {0, yMax}},
    AspectRatio -> resAspect,
    ImageSize   -> resImageSize
  ]
]

PlotResiduals[specs : {__Association}, opts : OptionsPattern[]] := Module[
  {th},
  th = ThemeFromOpts[{opts}];
  GraphicsGrid[
    {{PlotResidualMag[specs, Theme -> th], PlotResidualPhase[specs, Theme -> th]}},
    ImageSize  -> 900,
    Spacings   -> {0.5, 0.5},
    Background -> th["Bg"]
  ]
]

End[]
EndPackage[]
