(* ::Package:: *)

(* DCTPlots`PlotImpedance
   Nyquist and Bode plots for raw and fitted impedance.
   All functions take a single spec Association (or list) and return Graphics.

   Public API:
     PlotNyquist[spec]          Nyquist: -Im[Z] vs Re[Z], data + fit
     PlotBodeMag[spec]          Bode magnitude: |Z| vs f
     PlotBodePhase[spec]        Bode phase: angle(Z) vs f (degrees)
     PlotBode[spec]             GraphicsGrid: magnitude + phase side by side

   Theme-aware: default DCTPlots`$DCTTheme ("Dark" | "Publication"); override
   one call with Theme -> "Publication". Log axes use ScalingFunctions so the
   theme outward-tick generators apply.
*)

BeginPackage["DCTPlots`"]

PlotNyquist::usage =
  "PlotNyquist[spec] returns a Nyquist plot (-Im[Z] vs Re[Z]) overlaying \
raw data (dots) and DCT fit (line)."

PlotBodeMag::usage =
  "PlotBodeMag[spec] returns a log-log Bode magnitude plot (|Z| vs frequency)."

PlotBodePhase::usage =
  "PlotBodePhase[spec] returns a semilog Bode phase plot (angle in degrees vs frequency)."

PlotBode::usage =
  "PlotBode[spec] returns a GraphicsGrid with magnitude and phase Bode panels."

Begin["`Private`"]

bodeImageSize = 420;
bodeAspect    = 0.72;
fitStyle      = Directive[Orange, AbsoluteThickness[2.]];      (* data colour, theme-independent *)
dataStyleFor[fg_] := Directive[fg, Opacity[0.7], PointSize[0.012]];

(* ------------------------------------------------------------------ *)

Options[PlotNyquist] = {Theme -> Automatic};

PlotNyquist[spec_Association, opts : OptionsPattern[]] := Module[
  {th, fg, plotOpts, zData, zFit, dataPts, fitPts, xr, yr},

  th       = ThemeFromOpts[{opts}];
  fg       = th["Fg"];
  plotOpts = FilterRules[{opts}, Options[ListLinePlot]];

  zData   = spec["ZData"];
  zFit    = spec["ZFit"];
  dataPts = Transpose[{Re[zData], -Im[zData]}];
  fitPts  = Transpose[{Re[zFit],  -Im[zFit]}];
  xr      = MinMax[dataPts[[All, 1]]];
  yr      = MinMax[dataPts[[All, 2]]];

  Show[
    ListPlot[dataPts, PlotTheme -> "Default",
      PlotStyle -> dataStyleFor[fg], PlotMarkers -> {Automatic, 6}],
    ListLinePlot[fitPts, PlotTheme -> "Default", PlotStyle -> fitStyle],
    plotOpts,
    Sequence @@ ThemeChrome[th, 14, 1.1],
    FrameTicks  -> {{ThemeLinTicks[yr[[1]], yr[[2]], fg], None},
                    {ThemeLinTicks[xr[[1]], xr[[2]], fg], None}},
    FrameTicksStyle -> Directive[fg, 12],
    FrameLabel  -> {
      Style["Re[Z] (\[CapitalOmega])", 14, fg],
      Style["-Im[Z] (\[CapitalOmega])", 14, fg]
    },
    AspectRatio -> 1,
    ImageSize   -> bodeImageSize,
    PlotRange   -> All
  ]
]

(* ------------------------------------------------------------------ *)

Options[PlotBodeMag] = {Theme -> Automatic};

PlotBodeMag[spec_Association, opts : OptionsPattern[]] := Module[
  {th, fg, plotOpts, freq, zData, zFit, dataPts, fitPts, fMin, fMax, zMin, zMax},

  th       = ThemeFromOpts[{opts}];
  fg       = th["Fg"];
  plotOpts = FilterRules[{opts}, Options[ListLinePlot]];

  freq    = spec["FreqHz"];
  zData   = spec["ZData"];
  zFit    = spec["ZFit"];
  dataPts = Transpose[{freq, Abs[zData]}];
  fitPts  = Transpose[{freq, Abs[zFit]}];
  {fMin, fMax} = MinMax[freq];
  {zMin, zMax} = MinMax[Join[Abs[zData], Abs[zFit]]];

  (* single ListLinePlot (not Show) so FrameTicks positions stay in real coords
     under ScalingFunctions: data series = markers, fit series = line *)
  ListLinePlot[
    {dataPts, fitPts},
    plotOpts,
    Joined           -> {False, True},
    PlotMarkers      -> {{Automatic, 5}, None},
    PlotStyle        -> {dataStyleFor[fg], fitStyle},
    ScalingFunctions -> {"Log10", "Log10"},
    Sequence @@ ThemeChrome[th, 14, 1.1],
    FrameTicks  -> {{ThemeLogTicks[zMin, zMax, fg], None},
                    {ThemeLogTicks[fMin, fMax, fg], None}},
    FrameTicksStyle -> Directive[fg, 12],
    FrameLabel  -> {
      Style["Frequency (Hz)", 14, fg],
      Style["|Z| (\[CapitalOmega])", 14, fg]
    },
    AspectRatio -> bodeAspect,
    ImageSize   -> bodeImageSize,
    PlotRange   -> All
  ]
]

(* ------------------------------------------------------------------ *)

Options[PlotBodePhase] = {Theme -> Automatic};

PlotBodePhase[spec_Association, opts : OptionsPattern[]] := Module[
  {th, fg, plotOpts, freq, zData, zFit, dataPts, fitPts, fMin, fMax, pMin, pMax},

  th       = ThemeFromOpts[{opts}];
  fg       = th["Fg"];
  plotOpts = FilterRules[{opts}, Options[ListLinePlot]];

  freq    = spec["FreqHz"];
  zData   = spec["ZData"];
  zFit    = spec["ZFit"];
  dataPts = Transpose[{freq, (180. / Pi) Arg /@ zData}];
  fitPts  = Transpose[{freq, (180. / Pi) Arg /@ zFit}];
  {fMin, fMax} = MinMax[freq];
  {pMin, pMax} = MinMax[Join[dataPts[[All, 2]], fitPts[[All, 2]]]];

  ListLinePlot[
    {dataPts, fitPts},
    plotOpts,
    Joined           -> {False, True},
    PlotMarkers      -> {{Automatic, 5}, None},
    PlotStyle        -> {dataStyleFor[fg], fitStyle},
    ScalingFunctions -> {"Log10", None},
    Sequence @@ ThemeChrome[th, 14, 1.1],
    FrameTicks  -> {{ThemeLinTicks[pMin, pMax, fg], None},
                    {ThemeLogTicks[fMin, fMax, fg], None}},
    FrameTicksStyle -> Directive[fg, 12],
    FrameLabel  -> {
      Style["Frequency (Hz)", 14, fg],
      Style["Phase (\[Degree])", 14, fg]
    },
    AspectRatio -> bodeAspect,
    ImageSize   -> bodeImageSize,
    PlotRange   -> All
  ]
]

(* ------------------------------------------------------------------ *)

Options[PlotBode] = {Theme -> Automatic};

PlotBode[spec_Association, opts : OptionsPattern[]] := Module[
  {th},
  th = ThemeFromOpts[{opts}];
  GraphicsGrid[
    {{PlotBodeMag[spec, Theme -> th], PlotBodePhase[spec, Theme -> th]}},
    ImageSize  -> 900,
    Spacings   -> {0.5, 0.5},
    Background -> th["Bg"]
  ]
]

End[]
EndPackage[]
