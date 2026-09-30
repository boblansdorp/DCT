(* ::Package:: *)

(* DCTKKCheck`
   Lin-KK validation of EIS data in admittance space.
   Based on Schönleber et al., Electrochimica Acta 131, 20-27 (2014),
   adapted to use the Maxwell (RC-parallel) admittance basis consistent
   with the DCT ladder model.

   Model: Y_KK(omega) = Y_inf + Sum_mu g_mu * (i omega tau_mu) / (1 + i omega tau_mu)

   Re[Y_KK] = Y_inf + Sum_mu g_mu * (omega tau_mu)^2 / (1 + (omega tau_mu)^2)
   Im[Y_KK] =         Sum_mu g_mu * omega tau_mu    / (1 + (omega tau_mu)^2)

   Coefficients g_mu are unconstrained (unlike DCT which uses NNLS).
   Solved with LeastSquares on stacked [Re(Y); Im(Y)].

   Public API:
     KKCheck[data]               -> Association with residuals, ZFit, quality flags
     KKCheck[data, m]            -> same, m = number of Maxwell elements (default: auto)
     KKValidFreqRange[kkResult]  -> {fMin, fMax} largest window where |residuals| < threshold
     KKValidFreqRange[kkResult, threshold]
*)

BeginPackage["DCTKKCheck`"]

KKCheck::usage =
  "KKCheck[data] performs a Lin-KK consistency check in admittance space on an \
EIS data Association (output of ImportEIS). Returns Association with keys: \
FreqHz, Z, ZFit, YResRe, YResIm, RMSRe, RMSIm, M, OK."

KKValidFreqRange::usage =
  "KKValidFreqRange[kkResult, threshold:0.02] returns {fMin, fMax}, the \
largest contiguous frequency sub-range where both |YResRe| and |YResIm| stay \
below threshold."

Begin["`Private`"]

kkTaus[freqHz_List, m_Integer] :=
  10.^Subdivide[
    Log10[1. / (2. Pi Max[freqHz])],
    Log10[1. / (2. Pi Min[freqHz])],
    m - 1
  ]

(* Returns {Are, Aim}: N x (M+1) matrices for Re and Im of Y_KK *)
kkAdmittanceMatrix[omega_List, taus_List] := Module[
  {n, m, wt, are, aim},
  n = Length[omega];
  m = Length[taus];
  wt = Outer[Times, omega, taus];                       (* N x M: omega_i * tau_mu *)
  are = ArrayFlatten[{{ConstantArray[1., {n, 1}],
         wt^2 / (1. + wt^2)}}];                         (* N x (M+1) *)
  aim = ArrayFlatten[{{ConstantArray[0., {n, 1}],
         wt   / (1. + wt^2)}}];                         (* N x (M+1) *)
  {are, aim}
]

KKCheck[data_Association] :=
  KKCheck[data, Length[data["FreqHz"]]]

KKCheck[data_Association, m_Integer] := Module[
  {freqHz, z, y, omega, taus, are, aim, aFull, bFull, theta,
   yReFit, yImFit, yFit, zFit, yResRe, yResIm, rmsRe, rmsIm},

  freqHz = data["FreqHz"];
  z      = data["Z"];
  y      = 1. / z;
  omega  = 2. Pi freqHz;

  taus         = kkTaus[freqHz, m];
  {are, aim}   = kkAdmittanceMatrix[omega, taus];

  aFull  = Join[are, aim];
  bFull  = Join[Re[y], Im[y]];
  theta  = LeastSquares[aFull, bFull];

  yReFit = are . theta;
  yImFit = aim . theta;
  yFit   = yReFit + I yImFit;
  zFit   = 1. / yFit;

  yResRe = (Re[y] - yReFit) / Abs[y];
  yResIm = (Im[y] - yImFit) / Abs[y];
  rmsRe  = Sqrt[Mean[yResRe^2]];
  rmsIm  = Sqrt[Mean[yResIm^2]];

  <|"FreqHz" -> freqHz,
    "Z"      -> z,
    "ZFit"   -> zFit,
    "YResRe" -> yResRe,
    "YResIm" -> yResIm,
    "RMSRe"  -> rmsRe,
    "RMSIm"  -> rmsIm,
    "M"      -> m,
    "OK"     -> (Max[rmsRe, rmsIm] < 0.02)|>
]

KKValidFreqRange[kkResult_Association, threshold_ : 0.02] := Module[
  {freq, resRe, resIm, pass, trueRuns, bestRun},

  freq   = kkResult["FreqHz"];
  resRe  = kkResult["YResRe"];
  resIm  = kkResult["YResIm"];

  pass = MapThread[Abs[#1] < threshold && Abs[#2] < threshold &, {resRe, resIm}];

  trueRuns = Select[
    Split[Range[Length[pass]], pass[[#1]] === pass[[#2]] &],
    pass[[First[#]]] &
  ];

  If[trueRuns === {},
    {Min[freq], Max[freq]},
    bestRun = MaximalBy[trueRuns, Length][[1]];
    {freq[[First[bestRun]]], freq[[Last[bestRun]]]}
  ]
]

End[]
EndPackage[]
