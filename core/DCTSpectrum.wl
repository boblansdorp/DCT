(* ::Package:: *)

(* DCTSpectrum
   Top-level entry point: loads an EIS file and returns the best-fit DCT result.
   Outer Brent loop optimises Rs; inner NNLS solve is delegated to DCTKernel.

   Public API (exported into DCT` context for backward compatibility):
     DCTSpectrum[file]
     DCTSpectrum[file, opts...]

   Options:
     "BinsPerDecade"   -> 25         log-spaced tau bins per decade
     "FMinUse"         -> 0.5        lower frequency bound for fit [Hz]
     "FMaxUse"         -> 5000       upper frequency bound for fit [Hz]
     "TauMinFactor"    -> 0.1        tau grid lower extension factor
     "TauMaxFactor"    -> 10.        tau grid upper extension factor
     "LambdaND"        -> 0.01       dimensionless Tikhonov regularisation
     "WeightPower"     -> 0.75       admittance weight exponent (|Y|^p)
     "Debug"           -> False

   Returns Association with keys:
     Tau, g, C0, Rs, FreqHz, ZData, ZFit, YIntData, YIntFit,
     FinishTimeS, ObjZ, ObjY, RsFitFreqHz, RsFitZData
*)

BeginPackage["DCT`", {"DCTDataImport`", "DCTKernel`"}]

DCTSpectrum::usage =
  "DCTSpectrum[file] fits an EIS .txt file using the Maxwell admittance \
ladder (DCT method). Returns an Association with Tau, g, C0, Rs, and \
reconstructed impedance. See Options[DCTSpectrum] for tuning parameters."

Begin["`Private`"]

Options[DCTSpectrum] = {
  "BinsPerDecade" -> 25,
  "FMinUse"       -> 0.5,
  "FMaxUse"       -> 5000,
  "TauMinFactor"  -> 0.1,
  "TauMaxFactor"  -> 10.,
  "LambdaND"      -> 0.01,
  "WeightPower"   -> 0.75,
  "Debug"         -> False
}

DCTSpectrum[file_String, opts : OptionsPattern[]] := Catch @ Module[
  {
    debug, bpd, fMin, fMax, tauMinF, tauMaxF, lambdaND, wPow,
    dat, freqAll, zAll, timeS, finishTime,
    keepInner, freqInner, zInner,
    keepOuter, freqOuter, zOuter,
    rs0, rsLo, rsHi, rsWindow,
    objective, callCache, rsBest, refineSol, bestSolve,
    out
  },

  debug    = OptionValue["Debug"];
  bpd      = OptionValue["BinsPerDecade"];
  fMin     = N @ OptionValue["FMinUse"];
  fMax     = N @ OptionValue["FMaxUse"];
  tauMinF  = OptionValue["TauMinFactor"];
  tauMaxF  = OptionValue["TauMaxFactor"];
  lambdaND = OptionValue["LambdaND"];
  wPow     = OptionValue["WeightPower"];

  (* ---- Import ---- *)
  dat = DCTDataImport`ImportEIS[file, debug];
  If[FailureQ[dat],
    Throw[Failure["DCTSpectrum", <|"Message" -> "Import failed: " <> dat["Message"]|>]]
  ];

  freqAll  = dat["FreqHz"];
  zAll     = Developer`ToPackedArray[dat["Z"]];
  timeS    = dat["TimeS"];
  finishTime = If[AllTrue[timeS, MissingQ], Missing["NoTime"], Max[Select[timeS, NumericQ]]];

  If[debug, Print["[DCT] freq range: ", {Min[freqAll], Max[freqAll]}]];

  (* ---- Frequency windows ---- *)
  keepInner = (fMin <= # <= fMax) & /@ freqAll;
  freqInner = Pick[freqAll, keepInner];
  zInner    = Pick[zAll,    keepInner];

  (* Rs objective uses the same window (can differ in future extensions) *)
  keepOuter = keepInner;
  freqOuter = freqInner;
  zOuter    = zInner;

  If[Length[freqInner] < 5,
    Throw[Failure["DCTSpectrum",
      <|"Message" -> "Too few points in fit range [" <> ToString[fMin] <> ", " <> ToString[fMax] <> "] Hz"|>]]
  ];

  (* ---- Rs bounds for Brent ---- *)
  rs0 = DCTKernel`EstimateRs[freqAll, zAll, fMin, fMax];
  If[!NumericQ[rs0], rs0 = 0.];

  rsWindow = Select[Re[zOuter], NumericQ];
  rsLo = Max[0., Min[rsWindow] - 0.5 (Max[rsWindow] - Min[rsWindow])];
  rsHi = Max[rsLo + 1., Max[rsWindow] + 0.5 (Max[rsWindow] - Min[rsWindow])];

  If[debug, Print["[DCT] Rs0=", rs0, "  bounds=[", rsLo, ", ", rsHi, "]"]];

  (* ---- Outer Brent loop ---- *)
  callCache = <||>;

  objective[rs_?NumericQ] := Module[{key, sol, obj},
    key = ToString @ NumberForm[N[rs], {16, 6}];
    If[KeyExistsQ[callCache, key], Return[callCache[key]]];

    sol = DCTKernel`SolveLadderGivenRs[
      rs, freqInner, zInner, freqOuter, zOuter,
      lambdaND, bpd, tauMinF, tauMaxF, wPow
    ];
    obj = If[TrueQ[sol["OK"]], N[sol["ObjZ"]], Infinity];
    callCache[key] = obj;
    If[debug, Print["[DCT] Rs=", N[rs, 6], "  Obj=", obj]];
    obj
  ];

  refineSol = Quiet @ Check[
    FindMinimum[
      {objective[r], rsLo <= r <= rsHi},
      {{r, rs0}},
      Method -> "Brent",
      MaxIterations -> 40,
      WorkingPrecision -> MachinePrecision
    ],
    $Failed
  ];

  rsBest = If[
    MatchQ[refineSol, {_?NumericQ, {__Rule}}],
    r /. refineSol[[2]],
    rs0
  ];
  If[!NumericQ[rsBest], rsBest = rs0];
  rsBest = N[rsBest];

  If[debug, Print["[DCT] best Rs = ", rsBest]];

  (* ---- Final solve at best Rs ---- *)
  bestSolve = DCTKernel`SolveLadderGivenRs[
    rsBest, freqInner, zInner, freqOuter, zOuter,
    lambdaND, bpd, tauMinF, tauMaxF, wPow
  ];

  If[!TrueQ[bestSolve["OK"]],
    Throw[Failure["DCTSpectrum",
      <|"Message" -> "Inner solve failed at Rs=" <> ToString[rsBest]|>]]
  ];

  <|
    "Tau"          -> bestSolve["Tau"],
    "g"            -> bestSolve["g"],
    "C0"           -> bestSolve["C0"],
    "Rs"           -> bestSolve["Rs"],
    "FreqHz"       -> bestSolve["FreqHz"],
    "ZData"        -> bestSolve["ZData"],
    "ZFit"         -> bestSolve["ZFit"],
    "YIntData"     -> bestSolve["YIntData"],
    "YIntFit"      -> bestSolve["YIntFit"],
    "FinishTimeS"  -> finishTime,
    "ObjZ"         -> bestSolve["ObjZ"],
    "ObjY"         -> bestSolve["ObjY"],
    "RsFitFreqHz"  -> freqOuter,
    "RsFitZData"   -> zOuter
  |>
]

End[]
EndPackage[]
