(* ::Package:: *)

(* End-to-end test of DCT`DCTSpectrum on noise-free synthetic EIS data.
   Tests the full pipeline: file I/O (ImportEIS) + Brent Rs optimisation + NNLS fit.

   Synthetic model:
     Z = Rs + 1 / (i*omega*C0 + g1*dLog * i*omega/(1 + i*omega*tau1))

   C0 must be large enough that omega_max * C0 >> sum(g*dLog/tau) so that
   Re[Z] approaches Rs at the highest measured frequency — this is what
   makes the Rs estimate (Min[Re[Z]]) accurate and the Brent loop fast.

   This test deliberately avoids duplicating TestDCTKernel: it does not
   test SolveLadderGivenRs directly, only that the file-IO and outer Brent
   loop produce correct end-to-end output.
*)

Needs["DCTDataImport`"]
Needs["DCTKernel`"]
Needs["DCT`"]

Begin["DCTTests`Inversion`Private`"]

$pass = 0; $fail = 0;

checkTrue[name_String, cond_] :=
  If[TrueQ[cond],
    ($pass++; Style["\[Checkmark] " <> name, Darker[Green]]),
    ($fail++; Style["\[Times] " <> name, Darker[Red]])
  ]

check[name_String, got_, expected_, tol_: 5*10^-2] :=
  If[NumericQ[got] && NumericQ[expected] &&
     Abs[got - expected] / Max[Abs[expected], 10^-30] < tol,
    ($pass++; Style["\[Checkmark] " <> name, Darker[Green]]),
    ($fail++;
     Row[{Style["\[Times] " <> name <> ": ", Darker[Red]],
          "got = ", got, "   expected = ", expected}])
  ]

(* ==== Synthetic data ==== *)

rsTrue  = 50.;
c0True  = 1.0*^-5;   (* 10 µF — large enough that omega_max*C0 >> g*dLog/tau *)
tau1    = 0.01;       (* s — single Maxwell peak *)
g1      = 1.0*^-4;   (* F/decade *)
dLog    = 1. / 15;   (* matches default binsPerDecade = 15 *)

freqHz = Developer`ToPackedArray @ (10.^Range[-1., 4., 0.10]);
omega  = 2 Pi freqHz;

yInt   = I*omega*c0True + g1*dLog * (I*omega) / (1 + I*omega*tau1);
zSynth = rsTrue + 1 / yInt;

(* Write as a Gamry-style TSV that ImportEIS can parse *)
rows    = MapThread[{#1, #2, #3} &, {freqHz, Re[zSynth], -Im[zSynth]}];
tmpFile = CreateTemporary[];
Export[tmpFile,
  StringJoin @ Riffle[
    Prepend[
      StringJoin[Riffle[ToString /@ #, "\t"]] & /@ rows,
      "Freq/Hz\tRe(Z)/Ohm\t-Im(Z)/Ohm"
    ],
    "\n"
  ],
  "Text"
];

(* ==== Verify file roundtrip ==== *)

imported = DCTDataImport`ImportEIS[tmpFile];
Print @ checkTrue["ImportEIS returns Association",  AssociationQ[imported]];
Print @ checkTrue["ImportEIS FreqHz length",        Length[imported["FreqHz"]] == Length[freqHz]];
Print @ checkTrue["ImportEIS Z is numeric complex", VectorQ[imported["Z"], NumericQ]];

(* Roundtrip: Re[Z] should match what we wrote *)
Print @ checkTrue["ImportEIS Re[Z] roundtrip (1% tol)",
  Max[Abs[(Re[imported["Z"]] - Re[zSynth]) / Re[zSynth]]] < 0.01
];

(* ==== Full DCTSpectrum fit ==== *)

spec = Check[
  DCT`DCTSpectrum[tmpFile,
    "FMinUse"       -> 0.5,
    "FMaxUse"       -> 2000.,
    "LambdaND"      -> 0.001,
    "WeightPower"   -> 0.5,
    "BinsPerDecade" -> 15
  ],
  $Failed
];

DeleteFile[tmpFile];

Print @ checkTrue["DCTSpectrum returns Association",  AssociationQ[spec]];
If[!AssociationQ[spec],
  Print[Style["Aborting remaining checks: DCTSpectrum failed.", Darker[Red], Bold]];
  $fail += 5;
  ,
  Print @ checkTrue["Spec g vector non-negative",  VectorQ[spec["g"], # >= -10^-12 &]];
  Print @ checkTrue["Spec ZFit is numeric vector", VectorQ[spec["ZFit"], NumericQ]];

  (* Rs recovery — good because C0 is large so rs0 ~ rsTrue *)
  Print @ check["Rs within 5% of truth", spec["Rs"], rsTrue, 0.05];

  (* Fit quality — noise-free data, expect very tight residual *)
  relRes = Norm[spec["ZFit"] - spec["ZData"]] / Norm[spec["ZData"]];
  Print @ checkTrue["Relative |Z| residual < 2%", relRes < 0.02];

  (* Peak location *)
  gFit   = spec["g"];
  tauFit = spec["Tau"];
  iPeak  = First @ Ordering[gFit, -1];
  Print @ checkTrue["Peak tau within 0.4 decades of tau1",
    Abs[Log10[tauFit[[iPeak]]] - Log10[tau1]] < 0.4
  ];

  (* C0 is only tested loosely — it can partially alias with nearby Maxwell bins *)
  Print @ checkTrue["C0 is positive", TrueQ[spec["C0"] > 0]];
];

Print[""];
Print[Style[
  "Inversion: " <> ToString[$pass] <> " passed, " <> ToString[$fail] <> " failed.",
  If[$fail == 0, Darker[Green], Darker[Red]], 14, Bold
]];

End[]
