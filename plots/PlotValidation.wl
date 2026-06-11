(* ::Package:: *)

(* DCTPlots`PlotValidation
   Synthetic-recovery validation figure. A known g(k) is forward-modelled to an
   EIS spectrum (DCTKernel`SimulateEIS) and re-inverted (DCT`DCTSpectrumData);
   this overlays recovered vs simulated g(k) above the simulated Bode panels.

   Public API:
     PlotCombValidation[kTrue, gTrue, recovered]

   Theme-aware: default DCTPlots`$DCTTheme ("Dark" | "Publication").
*)

BeginPackage["DCTPlots`"]

PlotCombValidation::usage =
  "PlotCombValidation[kTrue, gTrue, recovered] returns a 3-panel validation \
figure: top = simulated g(k) (dashed) vs recovered g(k) (solid) from a DCT \
re-inversion; bottom = the simulated Bode magnitude and phase. `recovered` is a \
DCTSpectrumData / DCTSpectrum result (uses keys Tau, g, FreqHz, ZData)."

Begin["`Private`"]

Options[PlotCombValidation] = {Theme -> Automatic};

PlotCombValidation[kTrue_List, gTrue_List, recovered_Association,
                   opts : OptionsPattern[]] := Module[
  {th, fg, gScale, kRec, gRec, kMin, kMax, gMax, simStyle, recStyle, gkPanel, spec},

  th     = ThemeFromOpts[{opts}];
  fg     = th["Fg"];
  gScale = 1.*^9;                       (* F -> nF *)
  kRec   = 1. / recovered["Tau"];
  gRec   = recovered["g"];
  kMin   = Min[kTrue];  kMax = Max[kTrue];
  gMax   = Max[Join[gTrue, gRec]] gScale * 1.3;   (* headroom so the legend clears the peaks *)

  simStyle = Directive[Gray, Dashing[{0.025, 0.02}], AbsoluteThickness[2.5]];
  recStyle = Directive[fg, AbsoluteThickness[2.]];

  (* top: simulated (dashed) vs recovered (solid) g(k) *)
  gkPanel = ListLinePlot[
    {Transpose[{kTrue, gTrue gScale}], Transpose[{kRec, gRec gScale}]},
    ScalingFunctions -> {"Log10", None},
    Joined           -> True,
    PlotStyle        -> {simStyle, recStyle},
    Sequence @@ ThemeChrome[th, 14, 1.2],
    FrameTicks       -> {{ThemeLinTicks[0., gMax, fg], None},
                         {ThemeLogTicks[kMin, kMax, fg], None}},
    FrameTicksStyle  -> Directive[fg, 12],
    FrameLabel       -> {
      Style[Row[{"Electron-transfer rate, k (", Superscript["s", -1], ")"}], 14, fg],
      Style["g(k) (nF/decade)", 14, fg]
    },
    (* pin x to the simulated k-range; the recovered grid extends far past it
       (near-zero tails) and would otherwise stretch the axis weirdly *)
    PlotRange        -> {{kMin, kMax}, {0, gMax}},
    PlotRangePadding -> {{Scaled[0.02], Scaled[0.02]}, {Scaled[0.02], 0}},
    PlotLegends      -> Placed[LineLegend[
        {Directive[Gray, Dashing[{0.025, 0.02}]], fg}, {"Simulated", "Recovered"},
        LabelStyle     -> Directive[fg, 12],
        Background     -> th["Bg"],
        LegendFunction -> (Framed[#, Background -> th["Bg"],
          FrameStyle -> Directive[fg, AbsoluteThickness[1.1]]] &)],
      {0.88, 0.88}],
    ImageSize        -> 880,
    AspectRatio      -> 0.42];

  (* bottom: simulated Bode magnitude + phase (reuse the themed PlotBode) *)
  spec = <|"FreqHz" -> recovered["FreqHz"],
           "ZData"  -> recovered["ZData"],
           "ZFit"   -> recovered["ZData"]|>;

  Column[{gkPanel, PlotBode[spec, Theme -> th]}, Alignment -> Center, Spacings -> 1.5]
]

End[]
EndPackage[]
