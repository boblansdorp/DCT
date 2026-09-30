(* ::Package:: *)

(* Structural tests for DCTPlots`ThemeLogTicks / ThemeLinTicks, self-contained
   theme-coloured outward-tick generators. Assertions only; the visual
   dark-vs-publication comparison lives in runTests.wb cells. *)

Needs["DCTPlots`"]

Begin["DCTTests`PlotTicks`Private`"]

$pass = 0; $fail = 0;

ok[name_String, cond_] :=
  If[TrueQ[cond],
    ($pass++; Print @ Style["\[Checkmark] " <> name, Darker[Green]]),
    ($fail++; Print @ Style["\[Times] " <> name, Darker[Red]])
  ]

(* outward tick: 3rd element {0, len} with len > 0 (CustomTicks uses 0 inside);
   major carries a Style[_, color] label, minor an empty "". *)
outwardQ[t_] := MatchQ[t, {_, _, {0 | 0., _?Positive}, ___}]
majorQ[t_]   := outwardQ[t] && t[[2]] =!= ""
minorQ[t_]   := outwardQ[t] && t[[2]] === ""
positions[ticks_, q_] := Cases[ticks, t_ /; q[t] :> t[[1]]]

(* ── ThemeLogTicks ─────────────────────────────────────────────────── *)

logT = DCTPlots`ThemeLogTicks[1., 1000.];

ok["ThemeLogTicks: list of well-formed outward ticks",
  ListQ[logT] && AllTrue[logT, outwardQ]];
ok["ThemeLogTicks: has both major and minor ticks",
  AnyTrue[logT, majorQ] && AnyTrue[logT, minorQ]];
ok["ThemeLogTicks: a major tick per decade (1, 10, 100, 1000)",
  Max[Abs[Sort[positions[logT, majorQ]] - {1., 10., 100., 1000.}]] < 10.^-6];
ok["ThemeLogTicks: major labels are Style[_, Black] by default",
  AllTrue[Cases[logT, t_ /; majorQ[t] :> t[[2]]], MatchQ[#, Style[_, Black]] &]];

logTw = DCTPlots`ThemeLogTicks[1., 1000., White];
ok["ThemeLogTicks: colour arg recolours labels (White)",
  AllTrue[Cases[logTw, t_ /; majorQ[t] :> t[[2]]], MatchQ[#, Style[_, White]] &]];

(* ── ThemeLinTicks ─────────────────────────────────────────────────── *)

linT = DCTPlots`ThemeLinTicks[0., 1.];

ok["ThemeLinTicks: list of well-formed outward ticks",
  ListQ[linT] && AllTrue[linT, outwardQ]];
ok["ThemeLinTicks: has both major and minor ticks",
  AnyTrue[linT, majorQ] && AnyTrue[linT, minorQ]];
ok["ThemeLinTicks: all ticks within [0, 1]",
  AllTrue[positions[linT, outwardQ], 0. <= # <= 1. &]];

(* ── Plot-function integration: outward ticks ──────────────────────── *)
(* PlotKKResiduals used to be a ListLogLinearPlot with inward ticks; it now
   builds on ScalingFunctions so the theme outward-tick generators apply.
   Regression guard: assert its frame ticks really come out outward. *)

frameAxisTicks[p_, axis_] := Module[{g, ft},
  g  = If[Head[p] === Legended, First[p], p];
  ft = FirstCase[g, HoldPattern[FrameTicks -> v_] :> v, $Failed, Infinity];
  If[ft === $Failed, {}, ft[[axis, 1]]]   (* axis: 1 = left (y), 2 = bottom (x) *)
];

kkF   = 10.^Range[-1., 4., 0.1];
kkSyn = <|"FreqHz" -> kkF,
          "YResRe" -> RandomReal[{-0.01, 0.01}, Length[kkF]],
          "YResIm" -> RandomReal[{-0.01, 0.01}, Length[kkF]]|>;
kkPlot = DCTPlots`PlotKKResiduals[kkSyn];

ok["PlotKKResiduals renders (Graphics/Legended)",
  MatchQ[Head[kkPlot], Graphics | Legended]];
ok["PlotKKResiduals x-axis (frequency) ticks all outward",
  With[{bt = frameAxisTicks[kkPlot, 2]}, bt =!= {} && AllTrue[bt, outwardQ]]];
ok["PlotKKResiduals y-axis (residual) ticks all outward",
  With[{lt = frameAxisTicks[kkPlot, 1]}, lt =!= {} && AllTrue[lt, outwardQ]]];

histPlot = DCTPlots`PlotKKFMinHistogram[RandomReal[{0.1, 30.}, 400]];

ok["PlotKKFMinHistogram renders (Graphics)",
  MatchQ[Head[histPlot], Graphics | Legended]];
ok["PlotKKFMinHistogram x-axis (fMin) ticks all outward",
  With[{bt = frameAxisTicks[histPlot, 2]}, bt =!= {} && AllTrue[bt, outwardQ]]];
ok["PlotKKFMinHistogram y-axis (count) ticks all outward",
  With[{lt = frameAxisTicks[histPlot, 1]}, lt =!= {} && AllTrue[lt, outwardQ]]];

(* PlotVariancePowerLaw was a Show + ListLogLogPlot (inward ticks, Log[] epilog
   coords); now ScalingFunctions with Log10 epilog coords. *)
vplF     = 10.^Range[-1., 4., 0.1];
vplStats = <|"FreqHz" -> vplF, "YMean" -> 1.*^-3 vplF^0.5,
             "YVar" -> 1.*^-13 vplF, "PowerLawAlpha" -> 2.0, "PowerLawC" -> 1.*^-7|>;
vplPlot  = DCTPlots`PlotVariancePowerLaw[vplStats];

ok["PlotVariancePowerLaw renders",
  MatchQ[Head[vplPlot], Graphics | Legended]];
ok["PlotVariancePowerLaw x-axis ticks all outward",
  With[{bt = frameAxisTicks[vplPlot, 2]}, bt =!= {} && AllTrue[bt, outwardQ]]];
ok["PlotVariancePowerLaw y-axis ticks all outward",
  With[{lt = frameAxisTicks[vplPlot, 1]}, lt =!= {} && AllTrue[lt, outwardQ]]];

(* PlotLCurve: log x (roughness), linear y (narrow misfit range), both outward. *)
lcData = Table[{10.^-x, 5.9*^-4 + 0.02*^-4 x}, {x, 5., 8., 0.3}];
lcPlot = DCTPlots`PlotLCurve[lcData];

ok["PlotLCurve renders",
  MatchQ[Head[lcPlot], Graphics | Legended]];
ok["PlotLCurve x-axis (roughness, log) ticks all outward",
  With[{bt = frameAxisTicks[lcPlot, 2]}, bt =!= {} && AllTrue[bt, outwardQ]]];
ok["PlotLCurve y-axis (misfit, linear) ticks all outward",
  With[{lt = frameAxisTicks[lcPlot, 1]}, lt =!= {} && AllTrue[lt, outwardQ]]];

(* ── Converted plot functions: outward ticks (Sections 2-5) ────────── *)
(* These plots were ListLogLogPlot/ListLogLinearPlot/ListLogPlot (inward ticks);
   now ScalingFunctions so the theme outward generators apply on both axes. *)

outwardAxis[p_, axis_] :=
  With[{t = frameAxisTicks[p, axis]}, t =!= {} && AllTrue[t, outwardQ]];

synF = 10.^Range[-1., 4., 0.2];
synSpec[rs_, peakK_] := Module[{w = 2 Pi synF, z, tau, kk, gg},
  z   = rs + 500./(1 + I w 500. 1.*^-5);
  tau = 10.^Range[-4., 0., 0.1]; kk = 1./tau;
  gg  = Exp[-((Log10[kk] - Log10[peakK])^2)/(2 0.3^2)];
  <|"FreqHz" -> synF, "ZData" -> z (1 + 0.01 Sin[w]), "ZFit" -> z,
    "YIntData" -> (1/z) (1 + 0.005 Cos[w]), "YIntFit" -> 1/z,
    "Rs" -> rs, "C0" -> 1.*^-5, "Tau" -> tau, "g" -> gg|>];
synGoodSpecs = Flatten @ Table[
  <|"File" -> "C:/d/E" <> ToString[e] <> "/E" <> ToString[e] <> "_00" <> ToString[i] <> ".txt",
    "Spec" -> Append[synSpec[50. + i, 50. + 20 i], "FinishTimeS" -> i*3600.]|>,
  {e, 2}, {i, 3}];
synLSweep = Table[
  <|"File" -> "\[Lambda]=" <> ToString[N[10.^-x]], "Spec" -> synSpec[50., 100.]|>, {x, 1, 4}];
synSp1 = synGoodSpecs[[1, "Spec"]];

ok["PlotResidualMag outward (both axes)",
  outwardAxis[DCTPlots`PlotResidualMag[synSp1], 1] &&
  outwardAxis[DCTPlots`PlotResidualMag[synSp1], 2]];
ok["PlotResidualPhase outward (both axes)",
  outwardAxis[DCTPlots`PlotResidualPhase[synSp1], 1] &&
  outwardAxis[DCTPlots`PlotResidualPhase[synSp1], 2]];
ok["PlotBodeMag outward (both axes)",
  outwardAxis[DCTPlots`PlotBodeMag[synSp1], 1] &&
  outwardAxis[DCTPlots`PlotBodeMag[synSp1], 2]];
ok["PlotBodePhase outward (both axes)",
  outwardAxis[DCTPlots`PlotBodePhase[synSp1], 1] &&
  outwardAxis[DCTPlots`PlotBodePhase[synSp1], 2]];
ok["PlotNyquist outward (both axes)",
  outwardAxis[DCTPlots`PlotNyquist[synSp1], 1] &&
  outwardAxis[DCTPlots`PlotNyquist[synSp1], 2]];
ok["PlotCV outward (both axes)",
  With[{cv = DCTPlots`PlotCV[<|"FreqHz" -> synF, "CV" -> 0.05 + 0.01 Sin[Log[synF]]|>]},
    outwardAxis[cv, 1] && outwardAxis[cv, 2]]];
ok["PlotRsVsTime outward (both axes)",
  outwardAxis[DCTPlots`PlotRsVsTime[synGoodSpecs], 1] &&
  outwardAxis[DCTPlots`PlotRsVsTime[synGoodSpecs], 2]];
ok["PlotKPeakVsTime outward (both axes)",
  outwardAxis[DCTPlots`PlotKPeakVsTime[synGoodSpecs], 1] &&
  outwardAxis[DCTPlots`PlotKPeakVsTime[synGoodSpecs], 2]];
ok["PlotLambdaResiduals outward (both axes)",
  outwardAxis[DCTPlots`PlotLambdaResiduals[synLSweep], 1] &&
  outwardAxis[DCTPlots`PlotLambdaResiduals[synLSweep], 2]];

(* ── Summary ───────────────────────────────────────────────────────── *)

Print[Style[
  "PlotTicks: " <> ToString[$pass] <> " passed, " <> ToString[$fail] <> " failed.",
  If[$fail == 0, Darker[Green], Darker[Red]], 14, Bold]];

End[]
