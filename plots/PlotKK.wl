(* ::Package:: *)

(* DCTPlots`PlotKK
   Visualisation of Lin-KK residuals from KKCheck[].

   Public API:
     PlotKKResiduals[kkResult]             residuals for one spectrum
     PlotKKResiduals[kkResult, threshold]  same, custom threshold bands
*)

BeginPackage["DCTPlots`"]

PlotKKResiduals::usage =
  "PlotKKResiduals[kkResult, threshold:0.02] plots the relative admittance \
residuals (real and imaginary) vs frequency from KKCheck[]. Dashed lines mark \
the pass/fail threshold."

Begin["`Private`"]

$kkDarkBg    = GrayLevel[0.12];
$kkFgStyle   = Directive[White, AbsoluteThickness[1.1]];
$kkLblStyle  = Directive[White, 13, FontFamily -> "Arial"];
$kkReStyle   = Directive[RGBColor[0.35, 0.75, 1.],  AbsoluteThickness[1.8]];
$kkImStyle   = Directive[RGBColor[1.,   0.65, 0.15], AbsoluteThickness[1.8]];
$kkBandStyle = Directive[White, Opacity[0.35], Dashing[{0.015, 0.01}], AbsoluteThickness[1.]];

PlotKKResiduals[kkResult_Association, threshold_ : 0.02, opts : OptionsPattern[]] := Module[
  {freq, resRe, resIm, ptsRe, ptsIm, fMin, fMax, yMax},

  freq   = kkResult["FreqHz"];
  resRe  = kkResult["YResRe"];
  resIm  = kkResult["YResIm"];
  ptsRe  = Transpose[{freq, resRe}];
  ptsIm  = Transpose[{freq, resIm}];
  {fMin, fMax} = MinMax[freq];
  yMax   = Max[1.5 threshold, 1.2 Max[Abs[resRe], Abs[resIm]]];

  Show[
    ListLogLinearPlot[
      {ptsRe, ptsIm},
      PlotStyle    -> {$kkReStyle, $kkImStyle},
      Joined       -> True,
      PlotLegends  -> Placed[LineLegend[
        {$kkReStyle, $kkImStyle},
        {"\[Delta]Re(Y)/|Y|", "\[Delta]Im(Y)/|Y|"},
        LegendMarkerSize -> 20,
        LabelStyle       -> Directive[White, 12],
        Background       -> $kkDarkBg,
        LegendFunction   -> (Framed[#, Background -> $kkDarkBg,
          FrameStyle -> $kkFgStyle] &)],
        {0.82, 0.88}]
    ],
    (* threshold band lines *)
    Graphics[{
      $kkBandStyle,
      Line[{{Log[fMin], threshold},  {Log[fMax], threshold}}],
      Line[{{Log[fMin], -threshold}, {Log[fMax], -threshold}}]
    }],
    Frame       -> True,
    Axes        -> False,
    Background  -> $kkDarkBg,
    FrameStyle  -> $kkFgStyle,
    LabelStyle  -> $kkLblStyle,
    FrameLabel  -> {
      Style["Frequency (Hz)", 13, White],
      Style["Relative admittance residual", 13, White]
    },
    PlotRange   -> {All, {-yMax, yMax}},
    GridLines   -> {None, {0}},
    GridLinesStyle -> Directive[White, Opacity[0.2], AbsoluteThickness[0.7]],
    AspectRatio -> 0.65,
    ImageSize   -> 620,
    opts
  ]
]

End[]
EndPackage[]
