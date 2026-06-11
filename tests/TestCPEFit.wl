(* ::Package:: *)

(* Tests for CPEFit`
   Covers CPEImpedance (forward model) and FitCPE (5-parameter recovery of the
   Rs + (Cdl || (Rct + CPE)) circuit) against synthetic spectra, noise-free and
   with 1% Gaussian noise.
*)

Needs["CPEFit`"]

Begin["DCTTests`CPEFit`Private`"]

$pass = 0; $fail = 0;

checkTrue[name_String, cond_] :=
  If[TrueQ[cond],
    ($pass++; Style["\[Checkmark] " <> name, Darker[Green]]),
    ($fail++; Style["\[Times] " <> name, Darker[Red]])
  ]

check[name_String, got_, expected_, tol_ : 1*10^-2] :=
  If[NumericQ[got] && NumericQ[expected] &&
     Abs[got - expected] / Max[Abs[expected], 10^-30] < tol,
    ($pass++; Style["\[Checkmark] " <> name, Darker[Green]]),
    ($fail++;
     Row[{Style["\[Times] " <> name <> ": ", Darker[Red]],
       "got=", got, "  expected=", expected}])
  ]

(* ==== forward model sanity ==== *)

pTrue = {157., 2.*^-6, 5000., 5.*^-7, 0.88};   (* Rs, Cdl, Rct, Q, n *)
freq  = 10.^Range[-0.4, 2.6, 0.1];             (* ~0.4 to 398 Hz *)
zTrue = CPEFit`CPEImpedance[pTrue, freq];

Print @ checkTrue["CPEImpedance same length as freq", Length[zTrue] === Length[freq]];
Print @ checkTrue["CPEImpedance is complex",
  VectorQ[zTrue, NumericQ] && Max[Abs[Im[zTrue]]] > 0];

(* high-frequency limit -> Rs (real); low-freq -> blocking (|Z| large) *)
zHi = First @ CPEFit`CPEImpedance[pTrue, {1.*^7}];
Print @ checkTrue["high-freq Re -> Rs", Abs[Re[zHi] - pTrue[[1]]] / pTrue[[1]] < 0.05];
Print @ checkTrue["high-freq |Im| small", Abs[Im[zHi]] < pTrue[[1]]];
zLo = First @ CPEFit`CPEImpedance[pTrue, {1.*^-3}];
Print @ checkTrue["low-freq blocking (|Z| >> Rs)", Abs[zLo] > 50 pTrue[[1]]];

(* ==== noise-free 5-parameter recovery ==== *)

fit0 = CPEFit`FitCPE[freq, zTrue];
Print @ checkTrue["FitCPE returns Association", AssociationQ[fit0]];
Print @ check["recover Rs",  fit0["Rs"],  pTrue[[1]], 1*^-2];
Print @ check["recover Cdl", fit0["Cdl"], pTrue[[2]], 1*^-2];
Print @ check["recover Rct", fit0["Rct"], pTrue[[3]], 2*^-2];
Print @ check["recover Q",   fit0["Q"],   pTrue[[4]], 2*^-2];
Print @ check["recover n",   fit0["n"],   pTrue[[5]], 1*^-2];
Print @ checkTrue["noise-free RelResid < 0.1%", fit0["RelResid"] < 1.*^-3];
Print @ checkTrue["ZFit same length as data", Length[fit0["ZFit"]] === Length[freq]];

(* ==== recovery with 1% Gaussian noise ==== *)

SeedRandom[42];
zNoisy = zTrue (1. + RandomVariate[NormalDistribution[0., 0.01], Length[zTrue]]
                   + I RandomVariate[NormalDistribution[0., 0.01], Length[zTrue]]);
fitN = CPEFit`FitCPE[freq, zNoisy];
Print @ check["noisy recover Rs",  fitN["Rs"],  pTrue[[1]], 0.10];
Print @ check["noisy recover Rct", fitN["Rct"], pTrue[[3]], 0.20];
Print @ check["noisy recover n",   fitN["n"],   pTrue[[5]], 0.05];
Print @ checkTrue["noisy RelResid < 3%", fitN["RelResid"] < 0.03];

Print[""];
Print[Style["CPEFit: " <> ToString[$pass] <> " passed, " <> ToString[$fail] <> " failed.",
  If[$fail == 0, Darker[Green], Darker[Red]], 14, Bold]];

End[]
