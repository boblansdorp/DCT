(* ::Package:: *)

(* Tests for DCTKernel`
   Covers BuildTauGrid, SecondDifferenceMatrix, and SolveLadderGivenRs
   against synthetic data with a known single-RC response.
*)

Needs["NNLSFit`"]
Needs["DCTKernel`"]

Begin["DCTTests`Kernel`Private`"]

$pass = 0; $fail = 0;

checkTrue[name_String, cond_] :=
  If[TrueQ[cond],
    ($pass++; Style["\[Checkmark] " <> name, Darker[Green]]),
    ($fail++; Style["\[Times] " <> name, Darker[Red]])
  ]

approxEq[a_, b_, tol_ : 5*10^-3] :=
  NumericQ[a] && NumericQ[b] && Abs[a - b] / Max[Abs[b], 10^-30] < tol

check[name_String, got_, expected_, tol_ : 5*10^-3] :=
  If[approxEq[got, expected, tol],
    ($pass++; Style["\[Checkmark] " <> name, Darker[Green]]),
    ($fail++;
     Row[{Style["\[Times] " <> name <> ": ", Darker[Red]],
       "got=", got, "  expected=", expected}])
  ]

(* ==== BuildTauGrid ==== *)

taus = DCTKernel`BuildTauGrid[{0.1, 1., 10., 100., 1000.}, 10, 1., 1.];
Print @ checkTrue["TauGrid is packed real vector", Developer`PackedArrayQ[taus]];
Print @ checkTrue["TauGrid sorted ascending",      OrderedQ[taus]];
Print @ checkTrue["TauGrid min ~ 1/(2pi*1000)",
  Abs[Min[taus] - 1/(2 Pi 1000.)] / (1/(2 Pi 1000.)) < 0.5
];
Print @ checkTrue["TauGrid max ~ 1/(2pi*0.1)",
  Abs[Max[taus] - 1/(2 Pi 0.1)] / (1/(2 Pi 0.1)) < 0.5
];

(* ==== SecondDifferenceMatrix ==== *)

d2 = DCTKernel`SecondDifferenceMatrix[5];
Print @ checkTrue["D2 dims {3,5}",  Dimensions[d2] === {3, 5}];
Print @ checkTrue["D2 row sums zero", Max[Abs[Total[Normal[d2], {2}]]] < 10^-12];

(* ==== SolveLadderGivenRs: single Randle element synthetic test ==== *)
(*
   The Maxwell ladder represents Randle (series-RC) elements, not parallel RC.
   Each element: Y_k = (i*w) / (1 + i*w*tau_k), so Z_int = tau/(C) + 1/(i*w*C).
   Synthetic: Z = Rs + Rct + 1/(i*w*C)  [series Rs, charge-transfer Rct, cap C]
              Y_int = (i*w*C) / (1 + i*w*Rct*C)
   This is exactly one Maxwell ladder bin, so NNLS should find a sharp peak at
   tau_star = Rct * C.
*)

rsTrue  = 50.;
rctTrue = 1000.;
cTrue   = 1*10^-5;
tauTrue = rctTrue * cTrue;   (* 0.01 s  -- the DCT peak location *)

freqHz = 10.^Range[-1, 3, 0.2];  (* 0.1 Hz to 1 kHz *)
omega  = 2 Pi freqHz;

(* synthetic impedance: series Rs + Rct + C  (one Randle element) *)
zSynth = rsTrue + rctTrue + 1 / (I omega * cTrue);

(* fit range: exclude DC limit *)
fMin = 0.1; fMax = 1000.;
keep = (fMin <= # <= fMax) & /@ freqHz;
fInner = Pick[freqHz, keep];
zInner = Pick[zSynth, keep];

sol = DCTKernel`SolveLadderGivenRs[
  rsTrue, fInner, zInner, fInner, zInner,
  10^-3, 20, 1., 1., 0.75
];

Print @ checkTrue["SolveLadder OK flag",  TrueQ[sol["OK"]]];
Print @ checkTrue["SolveLadder g non-neg", VectorQ[sol["g"], # >= -10^-12 &]];
Print @ checkTrue["SolveLadder ZFit numeric", VectorQ[sol["ZFit"], NumericQ]];
Print @ checkTrue["SolveLadder same length as input",
  Length[sol["ZFit"]] == Length[fInner]
];

(* relative fit residual — regularisation introduces some smoothing so allow 20% *)
relRes = Norm[sol["ZFit"] - zInner] / Norm[zInner];
Print @ checkTrue["SolveLadder relative residual < 20%", relRes < 0.20];

(* the peak of g should be near tauTrue *)
gFit  = sol["g"];
tauFit = sol["Tau"];
iPeak = First[Ordering[gFit, -1]];
Print @ checkTrue["SolveLadder peak tau within half-decade of true",
  Abs[Log10[tauFit[[iPeak]]] - Log10[tauTrue]] < 0.6
];

Print[""];
Print[Style["DCTKernel: " <> ToString[$pass] <> " passed, " <> ToString[$fail] <> " failed.",
  If[$fail == 0, Darker[Green], Darker[Red]], 14, Bold]];

End[]
