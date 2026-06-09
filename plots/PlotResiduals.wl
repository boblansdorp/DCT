(* ::Package:: *)

(* DCTPlots`PlotResiduals
   Residual analysis plots comparing raw data to the DCT fit.

   Public API:
     PlotResidualMag[spec]     |(Z_data - Z_fit) / Z_data| vs frequency
     PlotResidualPhase[spec]   |phase_data - phase_fit| in degrees vs f
     PlotResiduals[spec]       GraphicsGrid: magnitude + phase side by side

   Default colour scheme is dark (white-on-dark) for wolfbook display.
   For light-background export pass: Background->White,
     FrameStyle->Directive[Black,AbsoluteThickness[1.1]],
     LabelStyle->Directive[Black,14,FontFamily->"Arial"]
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

$darkBg       = GrayLevel[0.12];
resFrameStyle = Directive[White, AbsoluteThickness[1.1]];
resLabelStyle = Directive[White, 14, FontFamily -> "Arial"];
resStyle      = Directive[Lighter[Blue, 0.4], AbsoluteThickness[1.8]];
resImageSize  = 420;
resAspect     = 0.72;

PlotResidualMag[spec_Association, opts : OptionsPattern[]] := Module[
  {freq, zData, zFit, pts},

  freq  = spec["FreqHz"];
  zData = spec["ZData"];
  zFit  = spec["ZFit"];
  pts   = Transpose[{freq, Abs[(zData - zFit) / zData]}];

  ListLogLogPlot[
    pts,
    Joined     -> True,
    PlotStyle  -> resStyle,
    Frame      -> True,
    Axes       -> False,
    Background -> $darkBg,
    FrameStyle -> resFrameStyle,
    LabelStyle -> resLabelStyle,
    FrameLabel -> {
      Style["Frequency (Hz)", 14, White],
      Style["|(Z_data - Z_fit)/Z_data|", 14, White]
    },
    PlotRange   -> {All, {10^-4, 1}},
    AspectRatio -> resAspect,
    ImageSize   -> resImageSize,
    opts
  ]
]

PlotResidualPhase[spec_Association, opts : OptionsPattern[]] := Module[
  {freq, zData, zFit, pts},

  freq  = spec["FreqHz"];
  zData = spec["ZData"];
  zFit  = spec["ZFit"];
  pts   = Transpose[{
    freq,
    Abs[(180. / Pi) Arg /@ zData - (180. / Pi) Arg /@ zFit]
  }];

  ListLogLinearPlot[
    pts,
    Joined     -> True,
    PlotStyle  -> resStyle,
    Frame      -> True,
    Axes       -> False,
    Background -> $darkBg,
    FrameStyle -> resFrameStyle,
    LabelStyle -> resLabelStyle,
    FrameLabel -> {
      Style["Frequency (Hz)", 14, White],
      Style["|phase residual| (\[Degree])", 14, White]
    },
    PlotRange   -> All,
    AspectRatio -> resAspect,
    ImageSize   -> resImageSize,
    opts
  ]
]

PlotResiduals[spec_Association, opts : OptionsPattern[]] :=
  GraphicsGrid[
    {{PlotResidualMag[spec], PlotResidualPhase[spec]}},
    ImageSize  -> 900,
    Spacings   -> {0.5, 0.5},
    Background -> $darkBg,
    opts
  ]

End[]
EndPackage[]
