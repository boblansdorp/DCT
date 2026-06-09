(* ::Package:: *)

(* Tests for DCTBiophysics
   Symbolic limit checks and parameter recovery on synthetic data.
*)

Needs["DCTBiophysics`"]

Begin["DCTTests`Biophysics`Private`"]

$pass = 0; $fail = 0;

approxEq[a_, b_, tol_ : 10^-4] :=
  VectorQ[{a, b}, NumericQ] && Abs[a - b] / Max[Abs[b], 10^-30] < tol

checkTrue[name_String, cond_] :=
  If[TrueQ[cond],
    ($pass++; Style["\[Checkmark] " <> name, Darker[Green]]),
    ($fail++; Style["\[Times] " <> name, Darker[Red]])
  ]

checkVal[name_String, got_, expected_, tol_ : 10^-4] :=
  checkTrue[name, approxEq[N[got], N[expected], tol]]

(* -- Symbolic form display -- *)

Print[Style["3-state model", Bold, 14]]
Print[Row[{Subscript["f","f"], " = ", $FF3State}]]
Print[""]
Print[Style["4-state model (with non-folding fraction)", Bold, 14]]
Print[Row[{Subscript["f","f"], " = ", $FF4State}]]
Print[""]

(* -- T1-T4: symbolic limit checks -- *)

ksTest = 0.2; kdTest = 50.; nfTest = 0.3;

Print @ checkVal["T1 3-state ff(C=0) = KS/(1+KS)",
  $FF3State /. {conc -> 0, KS -> ksTest, KD -> kdTest},
  ksTest / (1 + ksTest)]

Print @ checkVal["T2 3-state ff(C=large) -> 1",
  $FF3State /. {conc -> 10^8, KS -> ksTest, KD -> kdTest},
  1., 10^-3]

Print @ checkVal["T3 4-state ff(C=0) = (1-NF)*KS/(1+KS)",
  $FF4State /. {conc -> 0, KS -> ksTest, KD -> kdTest, NF -> nfTest},
  (1 - nfTest) * ksTest / (1 + ksTest)]

Print @ checkVal["T4 4-state ff(C=large) -> 1-NF",
  $FF4State /. {conc -> 10^8, KS -> ksTest, KD -> kdTest, NF -> nfTest},
  1. - nfTest, 10^-3]

(* -- T5-T8: parameter recovery on clean synthetic data -- *)

trueKS = 0.10; trueKD = 80.; trueNF = 0.35;
concPts = {1., 3., 10., 30., 100., 300., 1000.};

synthClean = {#, N[$FF4State /. {conc -> #, KS -> trueKS, KD -> trueKD, NF -> trueNF}]} & /@ concPts;

fit = DCTBiophysics`FitFoldedFraction[synthClean];

Print @ checkTrue["T5 fit returns Association", AssociationQ[fit]]
Print @ checkVal["T6 4-state KD recovery", fit["Model4", "KD"], trueKD, 0.05]
Print @ checkVal["T7 4-state NF recovery", fit["Model4", "NF"], trueNF, 0.05]
Print @ checkTrue["T8 3-state KD is larger (no NF model inflates KD)",
  fit["Model3", "KD"] > fit["Model4", "KD"]]

Print[""]
Print[Style["Biophysics: " <> ToString[$pass] <> " passed, " <> ToString[$fail] <> " failed.",
  If[$fail == 0, Darker[Green], Darker[Red]], 14, Bold]]

End[]
