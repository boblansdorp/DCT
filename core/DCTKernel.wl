(* ::Package:: *)

(* DCTKernel
   Pure mathematical kernel for DCT fitting: tau grid, Maxwell ladder matrix,
   Tikhonov regularization, and the inner NNLS solve at a fixed Rs.

   All functions are pure (no file I/O, no global state).

   Public API:
     BuildTauGrid[freqHz, binsPerDecade, tauMinFactor, tauMaxFactor]
       -> packed real vector of tau values

     EstimateRs[freqHz, z, fMin, fMax]
       -> scalar initial Rs guess

     SecondDifferenceMatrix[nBins]
       -> SparseArray, shape {nBins-2, nBins}

     SolveLadderGivenRs[rsCand, freqHzInner, zDataInner, freqHzOuter, zDataOuter, opts]
       -> Association with keys: OK, Rs, Tau, g, C0, ZFit, YIntData, YIntFit, ObjZ, ObjY
*)

BeginPackage["DCTKernel`", {"NNLSFit`"}]

BuildTauGrid::usage =
  "BuildTauGrid[freqHz, binsPerDecade, tauMinFactor, tauMaxFactor] returns a \
log-spaced tau grid spanning the frequency range of freqHz."

EstimateRs::usage =
  "EstimateRs[freqHz, z, fMin, fMax] returns Min[Re[Z]] over the frequency \
window [fMin, fMax] as an initial Rs guess."

SecondDifferenceMatrix::usage =
  "SecondDifferenceMatrix[n] returns the (n-2) x n second-difference matrix \
for Tikhonov smoothing."

SolveLadderGivenRs::usage =
  "SolveLadderGivenRs[rs, freqInner, zInner, freqOuter, zOuter, lambdaND, \
binsPerDecade, tauMinFactor, tauMaxFactor, weightPower] fits the Maxwell \
admittance ladder at fixed Rs and returns an Association with fit results."

Begin["`Private`"]

(* ------------------------------------------------------------------ *)
(* Tau grid                                                            *)
(* ------------------------------------------------------------------ *)

BuildTauGrid[freqHz_List, binsPerDecade_Integer,
             tauMinFactor_?NumericQ, tauMaxFactor_?NumericQ] := Module[
  {fMax, fMin, tauMin, tauMax, tauMinUse, tauMaxUse},

  fMax = Max[freqHz]; fMin = Min[freqHz];
  tauMin = 1 / (2 Pi fMax);
  tauMax = 1 / (2 Pi fMin);

  tauMinUse = If[tauMinFactor > 0, tauMin * tauMinFactor, tauMin];
  tauMaxUse = If[tauMaxFactor > 0, tauMax * tauMaxFactor, tauMax];
  If[tauMaxUse <= tauMinUse, tauMaxUse = tauMax];

  Developer`ToPackedArray @
    10.^Range[Log10[tauMinUse], Log10[tauMaxUse], 1. / binsPerDecade]
]

(* ------------------------------------------------------------------ *)
(* Rs initial guess                                                    *)
(* ------------------------------------------------------------------ *)

EstimateRs[freqHz_List, z_List, fMin_?NumericQ, fMax_?NumericQ] := Module[
  {keep, reVals},
  keep   = (fMin <= # <= fMax) & /@ freqHz;
  reVals = Select[Re[Pick[z, keep]], NumericQ];
  If[reVals === {}, 0., Max[0., Min[reVals]]]
]

(* ------------------------------------------------------------------ *)
(* Second-difference matrix for Tikhonov smoothing                    *)
(* ------------------------------------------------------------------ *)

SecondDifferenceMatrix[n_Integer] := SparseArray[
  Join[
    Table[{k, k}   ->  1, {k, 1, n - 2}],
    Table[{k, k+1} -> -2, {k, 1, n - 2}],
    Table[{k, k+2} ->  1, {k, 1, n - 2}]
  ],
  {n - 2, n}
]

(* ------------------------------------------------------------------ *)
(* Inner loop: solve Maxwell admittance ladder at fixed Rs             *)
(* ------------------------------------------------------------------ *)

SolveLadderGivenRs[
  rsCand_?NumericQ,
  freqHzInner_List, zDataInner_List,   (* DCT fit range *)
  freqHzOuter_List, zDataOuter_List,   (* Rs objective range *)
  lambdaND_?NumericQ,
  binsPerDecade_Integer,
  tauMinFactor_?NumericQ,
  tauMaxFactor_?NumericQ,
  weightPower_?NumericQ
] := Module[
  {
    omegaInner, omegaOuter,
    zInt, yInt, w,
    tau, nBins, dLog10,
    Kinner, Kouter, Kuse, aMat, aW, bW, aRI, bRI,
    d2, dFull, aCols, sA2, lambda, aAug, bAug, xBest,
    c0, gFit, dLogVec,
    yIntFitInner, zFitInner,
    KouterUse, yIntFitOuter,
    yDataOuter, wOuter,
    objZ, objY
  },

  omegaInner = Developer`ToPackedArray[2 Pi freqHzInner];
  omegaOuter = Developer`ToPackedArray[2 Pi freqHzOuter];

  zInt = zDataInner - rsCand;
  If[Min[Abs[zInt]] <= 10^-15,
    Return[<|"OK" -> False, "Message" -> "Z - Rs near zero", "ObjZ" -> Infinity|>]
  ];

  yInt   = Developer`ToPackedArray[1 / zInt];
  w = Developer`ToPackedArray[Abs[yInt]^weightPower];

  If[!VectorQ[w, NumericQ] ||
     !FreeQ[w, _ComplexInfinity | _DirectedInfinity | Indeterminate | Infinity],
    Return[<|"OK" -> False, "Message" -> "Weights non-finite", "ObjZ" -> Infinity|>]
  ];

  tau   = BuildTauGrid[freqHzInner, binsPerDecade, tauMinFactor, tauMaxFactor];
  nBins = Length[tau];
  dLog10 = Developer`ToPackedArray[ConstantArray[1. / binsPerDecade, nBins]];
  dLogVec = dLog10;

  (* Maxwell kernel matrix on inner range *)
  Kinner = Developer`ToPackedArray @ Table[
    (I omegaInner[[j]]) / (1 + I omegaInner[[j]] tau[[k]]),
    {j, Length[omegaInner]}, {k, nBins}
  ];
  Kuse = Kinner . SparseArray[DiagonalMatrix[dLogVec]];

  (* Design matrix: [iw*C0_column | Kuse] *)
  aMat = Join[Transpose[{I * omegaInner}], Kuse, 2];
  (* NNLS has no weight argument — it minimises ||A x - b||^2 with equal weights.
     To get weighted least squares Sum[w_i * r_i^2], pre-scale both A and b by
     Sqrt[w]: then NNLS squares Sqrt[w_i]*r_i and recovers w_i * r_i^2. *)
  aW   = DiagonalMatrix[Sqrt[w]] . aMat;
  bW   = Sqrt[w] * yInt;
  aRI  = Join[Re[aW], Im[aW]];
  bRI  = Join[Re[bW], Im[bW]];

  If[!MatrixQ[aRI, NumericQ] || !VectorQ[bRI, NumericQ],
    Return[<|"OK" -> False, "Message" -> "Design matrix non-numeric", "ObjZ" -> Infinity|>]
  ];

  (* Tikhonov regularization *)
  d2    = SecondDifferenceMatrix[nBins];
  dFull = ArrayFlatten[{{ConstantArray[0., {nBins - 2, 1}], d2}}];
  aCols = aRI[[All, 2 ;;]];
  sA2   = Median[Total[aCols^2, {1}]];
  lambda = lambdaND * sA2;

  aAug = Join[aRI, Sqrt[lambda] * N[Normal[dFull]], 1];
  bAug = Join[bRI, ConstantArray[0., nBins - 2]];

  If[!MatrixQ[aAug, NumericQ] || Dimensions[aAug][[1]] =!= Length[bAug] ||
     !FreeQ[aAug, _ComplexInfinity | _DirectedInfinity | Indeterminate | Infinity],
    Return[<|"OK" -> False, "Message" -> "Augmented system invalid", "ObjZ" -> Infinity|>]
  ];

  xBest = Quiet @ Check[NNLSFit`NNLS[aAug, bAug], $Failed];

  If[xBest === $Failed || !VectorQ[xBest, NumericQ] || Length[xBest] =!= nBins + 1,
    Return[<|"OK" -> False, "Message" -> "NNLS failed", "ObjZ" -> Infinity|>]
  ];

  c0   = xBest[[1]];
  gFit = xBest[[2 ;;]];

  (* Reconstruct fit on inner range *)
  yIntFitInner = (I omegaInner) * c0 + Kuse . gFit;
  If[Min[Abs[yIntFitInner]] <= 10^-30,
    Return[<|"OK" -> False, "Message" -> "Fit admittance too small", "ObjZ" -> Infinity|>]
  ];
  zFitInner = rsCand + 1 / yIntFitInner;

  (* Reconstruct on outer (Rs objective) range *)
  Kouter = Developer`ToPackedArray @ Table[
    (I omegaOuter[[j]]) / (1 + I omegaOuter[[j]] tau[[k]]),
    {j, Length[omegaOuter]}, {k, nBins}
  ];
  KouterUse = Kouter . SparseArray[DiagonalMatrix[dLogVec]];
  yIntFitOuter = (I omegaOuter) * c0 + KouterUse . gFit;

  If[Min[Abs[yIntFitOuter]] <= 10^-30,
    Return[<|"OK" -> False, "Message" -> "Outer fit admittance too small", "ObjZ" -> Infinity|>]
  ];

  yDataOuter = Developer`ToPackedArray[1 / (zDataOuter - rsCand)];
  If[!VectorQ[yDataOuter, NumericQ] ||
     !FreeQ[yDataOuter, _ComplexInfinity | _DirectedInfinity | Indeterminate | Infinity],
    Return[<|"OK" -> False, "Message" -> "Outer admittance data invalid", "ObjZ" -> Infinity|>]
  ];

  wOuter = Developer`ToPackedArray[Abs[yDataOuter]^weightPower];
  objZ = Total[wOuter * (Re[yDataOuter - yIntFitOuter]^2 + Im[yDataOuter - yIntFitOuter]^2)];
  objY = Total[Re[yInt - yIntFitInner]^2 + Im[yInt - yIntFitInner]^2];

  <|
    "OK"        -> True,
    "Message"   -> "",
    "Rs"        -> rsCand,
    "Tau"       -> tau,
    "g"         -> gFit,
    "C0"        -> c0,
    "FreqHz"    -> freqHzInner,
    "ZData"     -> zDataInner,
    "ZFit"      -> zFitInner,
    "YIntData"  -> yInt,
    "YIntFit"   -> yIntFitInner,
    "ObjZ"      -> objZ,
    "ObjY"      -> objY
  |>
]

End[]
EndPackage[]
