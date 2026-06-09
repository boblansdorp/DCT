(* ::Package:: *)

(* Tests for NNLSFit`NNLS
   Verifies correctness against known analytic solutions.
*)

Needs["NNLSFit`"]

Begin["DCTTests`NNLS`Private`"]

$pass = 0; $fail = 0;

approxEq[a_, b_, tol_ : 10^-6] :=
  VectorQ[{a, b}, NumericQ] && Abs[a - b] / Max[Abs[b], 10^-30] < tol

checkVec[name_String, got_, expected_, tol_ : 10^-5] := Module[{ok},
  ok = VectorQ[got, NumericQ] && Length[got] == Length[expected] &&
       And @@ MapThread[approxEq[#1, #2, tol] &, {got, expected}];
  If[ok,
    $pass++;
    Style["\[Checkmark] " <> name, Darker[Green]],
    $fail++;
    Row[{Style["\[Times] " <> name <> ": ", Darker[Red]], "got=", got, "  expected=", expected}]
  ]
]

checkTrue[name_String, cond_] :=
  If[TrueQ[cond],
    ($pass++; Style["\[Checkmark] " <> name, Darker[Green]]),
    ($fail++; Style["\[Times] " <> name, Darker[Red]])
  ]

(* T1: unconstrained solution is non-negative -> NNLS == least squares *)
(* A = identity, f = [3, 1, 2] -> x = [3, 1, 2] *)
A1 = IdentityMatrix[3];
f1 = {3., 1., 2.};
x1 = NNLSFit`NNLS[A1, f1];
Print @ checkVec["T1 identity: x=f when f>=0", x1, f1];

(* T2: one component of f is negative -> that component must be zeroed *)
(* A = identity, f = [1, -2, 3] -> x = [1, 0, 3] *)
A2 = IdentityMatrix[3];
f2 = {1., -2., 3.};
x2 = NNLSFit`NNLS[A2, f2];
Print @ checkVec["T2 negative f component -> zero", x2, {1., 0., 3.}];

(* T3: overdetermined system with known non-negative solution *)
(* A = {{1,0},{0,1},{1,1}}, f = {2,3,5} -> x = {2,3} exactly *)
A3 = {{1.,0.},{0.,1.},{1.,1.}};
f3 = {2., 3., 5.};
x3 = NNLSFit`NNLS[A3, f3];
Print @ checkVec["T3 overdetermined with exact solution", x3, {2., 3.}];

(* T4: result must be non-negative in all components *)
A4 = RandomReal[{0.1, 2.}, {10, 5}];
f4 = RandomReal[{-1., 1.}, 10];
x4 = NNLSFit`NNLS[A4, f4];
Print @ checkTrue["T4 all components >= 0", VectorQ[x4, # >= -10^-12 &]];
Print @ checkTrue["T4 returns real vector",  VectorQ[x4, NumericQ]];
Print @ checkTrue["T4 correct length",       Length[x4] == 5];

(* T5: residual at solution is at most as large as at zero *)
(* This verifies optimality: ||Ax - f|| <= ||f|| *)
A5 = RandomReal[{0.1, 1.}, {8, 4}];
f5 = RandomReal[{0.5, 2.}, 8];
x5 = NNLSFit`NNLS[A5, f5];
Print @ checkTrue["T5 residual <= ||f||",
  Norm[A5 . x5 - f5] <= Norm[f5] + 10^-8
];

Print[""];
Print[Style["NNLS: " <> ToString[$pass] <> " passed, " <> ToString[$fail] <> " failed.",
  If[$fail == 0, Darker[Green], Darker[Red]], 14, Bold]];

End[]
