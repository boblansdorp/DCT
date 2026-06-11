(* ::Package:: *)

(* NNLSFit
   Non-negative least squares: minimise ||A x - f||^2 subject to x >= 0.
   Uses the active-set / Lawson-Hanson algorithm with QR decomposition.

   Adapted from the Mathematica implementation by Michael Woodhams, posted to
   comp.soft-sys.math.mathematica (2 Oct 2003), of the Lawson & Hanson active-set
   algorithm ("Solving Least Squares Problems", Prentice-Hall 1974; SIAM 1995).
   Woodhams placed his code in the public domain and asks that his authorship be
   acknowledged.
     Original: https://groups.google.com/g/comp.soft-sys.math.mathematica/c/cHFiKQ8ssaI
     Repost:   https://mathematica.stackexchange.com/questions/269727

   Public API:
     NNLS[A, f]  ->  x  (real vector, length = Dimensions[A][[2]])
*)

BeginPackage["NNLSFit`"]

NNLS::usage =
  "NNLS[A, f] solves min_{x>=0} ||A.x - f||^2 using the Lawson-Hanson \
active-set algorithm. Returns a non-negative real vector x."

Begin["`Private`"]

activeIndices[zeroed_]   := Select[Range[Length[zeroed]], zeroed[[#]] == 1 &]
inactiveIndices[zeroed_] := Select[Range[Length[zeroed]], zeroed[[#]] == 0 &]

subMatrix[A_, cols_List] := Transpose[Transpose[A][[cols]]]

NNLS[A_?MatrixQ, f_?VectorQ] := Module[
  {n, x, zeroed, w, t, z, alpha, Ap, QR, compressedZ,
   positiveSet, zeroedSet, toBeZeroed},

  n = Dimensions[A][[2]];
  x      = ConstantArray[0., n];
  zeroed = ConstantArray[1, n];
  w      = Transpose[A] . (f - A . x);

  While[
    zeroedSet = activeIndices[zeroed];
    zeroedSet =!= {} && Max[w[[zeroedSet]]] > 0,

    (* move largest-gradient zeroed variable into positive set *)
    t = First[First[Position[w * zeroed, Max[w * zeroed]]]];
    zeroed[[t]] = 0;

    (* inner loop: ensure z >= 0 on positive set *)
    While[
      positiveSet  = inactiveIndices[zeroed];
      Ap           = subMatrix[A, positiveSet];
      QR           = QRDecomposition[Ap];
      compressedZ  = LinearSolve[QR[[2]], QR[[1]] . f];

      z            = ConstantArray[0., n];
      z[[positiveSet]] = compressedZ;

      Min[z] < 0,

      (* step toward feasibility *)
      alpha = Min[
        x[[#]] / (x[[#]] - z[[#]]) & /@
          Select[positiveSet, zeroed[[#]] == 0 && z[[#]] < 0 &]
      ];
      x = x + alpha * (z - x);

      toBeZeroed = Select[positiveSet, Abs[x[[#]]] < 10^-13 &];
      zeroed[[toBeZeroed]] = 1;
      x[[toBeZeroed]] = 0.
    ];

    x = z;
    w = Transpose[A] . (f - A . x)
  ];

  x
]

End[]
EndPackage[]
