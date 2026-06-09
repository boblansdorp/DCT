(* ::Package:: *)

(* DCTPlots`PlotImpedance
   Nyquist and Bode plots for raw and fitted impedance.
   All functions take a single spec Association (or list) and return Graphics.

   Public API:
     PlotNyquist[spec]          Nyquist: -Im[Z] vs Re[Z], data + fit
     PlotBodeMag[spec]          Bode magnitude: |Z| vs f
     PlotBodePhase[spec]        Bode phase: angle(Z) vs f (degrees)
     PlotBode[spec]             GraphicsGrid: magnitude + phase side by side

   Default colour scheme is dark (white-on-dark) for wolfbook display.
   For light-background export pass: Background->White,
     FrameStyle->Directive[Black,AbsoluteThickness[1.1]],
     LabelStyle->Directive[Black,14,FontFamily->"Arial"]
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

bodeLabelStyle = Directive[White, 14, FontFamily -> "Arial"];
bodeFrameStyle = Directive[White, AbsoluteThickness[1.1]];
bodeImageSize  = 420;
bodeAspect     = 0.72;
$darkBg        = GrayLevel[0.12];
dataStyle      = Directive[White, Opacity[0.7], PointSize[0.012]];
fitStyle       = Directive[Orange, AbsoluteThickness[2.]];

(* ------------------------------------------------------------------ *)

PlotNyquist[spec_Association, opts : OptionsPattern[]] := Module[
  {zData, zFit, dataPts, fitPts},

  zData   = spec["ZData"];
  zFit    = spec["ZFit"];
  dataPts = Transpose[{Re[zData], -Im[zData]}];
  fitPts  = Transpose[{Re[zFit],  -Im[zFit]}];

  Show[
    ListPlot[dataPts,
      PlotStyle   -> dataStyle,
      PlotMarkers -> {Automatic, 6}
    ],
    ListLinePlot[fitPts, PlotStyle -> fitStyle],
    Frame      -> True,
    Axes       -> False,
    Background -> $darkBg,
    FrameStyle -> bodeFrameStyle,
    LabelStyle -> bodeLabelStyle,
    FrameLabel -> {
      Style["Re[Z] (\[CapitalOmega])", 14, White],
      Style["-Im[Z] (\[CapitalOmega])", 14, White]
    },
    AspectRatio -> 1,
    ImageSize   -> bodeImageSize,
    PlotRange   -> All,
    opts
  ]
]

(* ------------------------------------------------------------------ *)

PlotBodeMag[spec_Association, opts : OptionsPattern[]] := Module[
  {freq, zData, zFit, dataPts, fitPts},

  freq    = spec["FreqHz"];
  zData   = spec["ZData"];
  zFit    = spec["ZFit"];
  dataPts = Transpose[{freq, Abs[zData]}];
  fitPts  = Transpose[{freq, Abs[zFit]}];

  Show[
    ListLogLogPlot[dataPts, PlotStyle -> dataStyle, Joined -> False],
    ListLogLogPlot[fitPts,  PlotStyle -> fitStyle,  Joined -> True],
    Frame      -> True,
    Axes       -> False,
    Background -> $darkBg,
    FrameStyle -> bodeFrameStyle,
    LabelStyle -> bodeLabelStyle,
    FrameLabel -> {
      Style["Frequency (Hz)", 14, White],
      Style["|Z| (\[CapitalOmega])", 14, White]
    },
    AspectRatio -> bodeAspect,
    ImageSize   -> bodeImageSize,
    PlotRange   -> All,
    opts
  ]
]

(* ------------------------------------------------------------------ *)

PlotBodePhase[spec_Association, opts : OptionsPattern[]] := Module[
  {freq, zData, zFit, dataPts, fitPts},

  freq    = spec["FreqHz"];
  zData   = spec["ZData"];
  zFit    = spec["ZFit"];
  dataPts = Transpose[{freq, (180. / Pi) Arg /@ zData}];
  fitPts  = Transpose[{freq, (180. / Pi) Arg /@ zFit}];

  Show[
    ListLogLinearPlot[dataPts, PlotStyle -> dataStyle, Joined -> False],
    ListLogLinearPlot[fitPts,  PlotStyle -> fitStyle,  Joined -> True],
    Frame      -> True,
    Axes       -> False,
    Background -> $darkBg,
    FrameStyle -> bodeFrameStyle,
    LabelStyle -> bodeLabelStyle,
    FrameLabel -> {
      Style["Frequency (Hz)", 14, White],
      Style["Phase (\[Degree])", 14, White]
    },
    AspectRatio -> bodeAspect,
    ImageSize   -> bodeImageSize,
    PlotRange   -> All,
    opts
  ]
]

(* ------------------------------------------------------------------ *)

PlotBode[spec_Association, opts : OptionsPattern[]] :=
  GraphicsGrid[
    {{PlotBodeMag[spec], PlotBodePhase[spec]}},
    ImageSize  -> 900,
    Spacings   -> {0.5, 0.5},
    Background -> $darkBg,
    opts
  ]

End[]
EndPackage[]
