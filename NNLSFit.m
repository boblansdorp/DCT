(* ::Package:: *)

BeginPackage["NNLSFit`"]

NNLS::usage = "NNLS[A, f] solves the non-negative least squares problem minimizing ||Ax - f||^2 with x >= 0.";

Begin["`Private`"]

bitsToIndices[v_] := Select[Range[Length[v]], v[[#]] == 1 &];

NNLS[A_, f_] := Module[
  {x, zeroed, w, t, Ap, z, q, \[Alpha], zeroedSet, positiveSet,
   toBeZeroed, compressedZ, Q, R},

  zeroedSet := bitsToIndices[zeroed];
  positiveSet := bitsToIndices[1 - zeroed];

  x = ConstantArray[0., Length@A[[1]]];
  zeroed = ConstantArray[1, Length@x];
  w = Transpose[A] . (f - A . x);

  While[zeroedSet =!= {} && Max[w[[zeroedSet]]] > 0,
    t = First@First@Position[w*zeroed, Max[w*zeroed]];
    zeroed[[t]] = 0;

    Ap = Transpose[Transpose[A][[positiveSet]]];
    {Q, R} = QRDecomposition[Ap];
    compressedZ = Inverse[R] . Q . f;

    z = ConstantArray[0., Length@x];
    z[[positiveSet]] = compressedZ;

    While[Min[z] < 0,
      \[Alpha] = Infinity;
      Do[
        If[zeroed[[q]] == 0 && z[[q]] < 0,
          \[Alpha] = Min[\[Alpha], x[[q]]/(x[[q]] - z[[q]])];
        ],
        {q, Length@x}
      ];

      x = x + \[Alpha] (z - x);
      toBeZeroed = Select[positiveSet, Abs[x[[#]]] < 10^-13 &];
      zeroed[[toBeZeroed]] = 1;
      x[[toBeZeroed]] = 0;

      Ap = Transpose[Transpose[A][[positiveSet]]];
      {Q, R} = QRDecomposition[Ap];
      compressedZ = Inverse[R] . Q . f;
      z = ConstantArray[0., Length@x];
      z[[positiveSet]] = compressedZ;
    ];

    x = z;
    w = Transpose[A] . (f - A . x);
  ];

  x
]

End[]

Print[Style["\:2714 NNLSFit package loaded successfully.", Darker[Green]]];

EndPackage[]

