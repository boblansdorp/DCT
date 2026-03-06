(* ::Package:: *)

(*  DCT.wl
    Bare-bones Wolfram Language package to compute a Distribution of Capacitance Timescales (DCT)
    spectrum from a Gamry/typical EIS .txt export that is semicolon-separated.

    Expected data columns (minimum):
    "Frequency (Hz)"
     "Z' (\[CapitalOmega])"
     "-Z'' (\[CapitalOmega])"

    Primary entry point:
      DCT`DCTSpectrum[file]

    Returns an Association containing:
      "Tau"   -> tau grid [s]
      "g"     -> capacitance density [F/decade] (default convention)
      "C0"    -> parallel capacitance term [F] (handles tau -> 0 / very fast processes)
      "Rs"    -> estimated series resistance [Ohm]
      "ZData" -> complex impedance data
      "ZFit"  -> reconstructed best-fit impedance
      "FreqHz"-> frequency vector [Hz]
*)
(* Load NNLSFit package from the same folder as this DCT .wl file *)
Module[{here, nnlsPathWL, nnlsPathM},
  here = DirectoryName[$InputFileName];
  nnlsPathM = "C:\\Users\\Bob Lansdorp\\Documents\\DCT\\NNLSFit.m";

  If[FileExistsQ[nnlsPathWL],
    Get[nnlsPathWL],
    If[FileExistsQ[nnlsPathM],
      Get[nnlsPathM],
      Print["NNLSFit not found next to DCT file. Expected: ", nnlsPathWL, " or ", nnlsPathM];
      Abort[];
    ]
  ];
];


BeginPackage["DCT`"];

DCTSpectrum::usage =
"DCTSpectrum[file] imports an EIS .txt file and returns an Association with the best-fit DCT spectrum.\n" <>
"Options let you adjust tau grid density, regularization, and Rs estimation.";

Begin["`Private`"];
(* debugging *)
dbg[flag_, msg_, expr_] := If[TrueQ[flag], Print[msg, ": ", expr]];

assert[cond_, msg_] := If[Not@TrueQ[cond], Throw[Failure["DCTDebug", <|"Message" -> msg|>]]];

dbg[flag_, label_, expr_] := If[TrueQ[flag], Print["[DCT DEBUG] ", label, " = ", expr]];
assert[cond_, msg_] := If[Not@TrueQ[cond], Throw[Failure["DCTDebug", <|"Message" -> msg|>]]];

(* end debug *)

(* ========================= *)
(* Helpers                   *)
(* ========================= *)

toNumber[s_String] := Quiet@Check[ToExpression[StringTrim[s]], Missing["Bad"]];
toNumber[x_] := x;

(* Import and parse semicolon-separated EIS text file.
   Minimal assumptions: at least 4 semicolon-separated fields per row. *)
importEISTxt[file_String] := Module[
 {rawLines, lines, header, hasHeader, splitHeader, colMap,
  freqCol, zreCol, zimCol, timeCol,
  splitRows, numRows, good,
  freqHz, zre, zimNeg, z, timeS},

 rawLines = Import[file, "Lines"];
 rawLines = Select[rawLines, StringTrim[#] =!= "" &];

 header = First[rawLines];
 hasHeader = StringContainsQ[header, "Frequency"];

 lines = If[hasHeader, Rest[rawLines], rawLines];

 If[hasHeader,
  
  splitHeader = StringTrim /@ StringSplit[header, ";"];
  colMap = AssociationThread[splitHeader -> Range[Length[splitHeader]]];

  freqCol = Lookup[colMap, "Frequency (Hz)", Missing["NoCol"]];
  zreCol  = Lookup[colMap, "Z' (\[CapitalOmega])", Missing["NoCol"]];
  zimCol  = Lookup[colMap, "-Z'' (\[CapitalOmega])", Missing["NoCol"]];
  timeCol = Lookup[colMap, "Time (s)", Missing["NoCol"]];

  If[MemberQ[{freqCol, zreCol, zimCol}, Missing["NoCol"]],
   Throw[
    Failure["DCTImport",
     <|"Message" -> "Header found but required columns not located."|>
    ]
   ]
  ],

  (* legacy format *)
  freqCol = 2; 
  zreCol = 3; 
  zimCol = 4;
  timeCol = Missing["NoCol"];
 ];

 splitRows = StringSplit[#, ";"] & /@ lines;
 numRows = (toNumber /@ #) & /@ splitRows;

 good = Select[
   numRows,
   (Length[#] >= Max[freqCol, zreCol, zimCol] &&
      NumericQ[#[[freqCol]]] &&
      NumericQ[#[[zreCol]]] &&
      NumericQ[#[[zimCol]]]) &
   ];

 freqHz = good[[All, freqCol]];
 zre    = good[[All, zreCol]];
 zimNeg = good[[All, zimCol]];

 z = zre + I*(-zimNeg);

 timeS =
  If[timeCol === Missing["NoCol"],
   ConstantArray[Missing["NoTime"], Length[freqHz]],
   good[[All, timeCol]]
   ];

 With[{ord = Ordering[freqHz]},
  <|
   "FreqHz" -> freqHz[[ord]],
   "Z" -> z[[ord]],
   "TimeS" -> timeS[[ord]]
   |>
 ]
];


(* Estimate Rs from points in a specified frequency window.
   Returns the median Re[Z] over minFreqRsFit <= f <= maxFreqRsFit. *)
estimateRs[freqHz_List, z_List, minFreqRsFit_?NumericQ, maxFreqRsFit_?NumericQ] := Module[
  {keep, top, rs},
  keep = (minFreqRsFit <= # <= maxFreqRsFit) & /@ freqHz;
  top = Re[Pick[z, keep]];
  top = Select[top, NumericQ];
  If[top === {}, Return[0.]];
  rs = Median[top];
  If[NumericQ[rs] && rs >= 0., rs, 0.]
];
(* Build a log-spaced tau grid. *)
buildTauGrid[freqHz_List, binsPerDecade_Integer, tauMinFactor_?NumericQ, tauMaxFactor_?NumericQ] := Module[
	{fMax, fMin, tauMin, tauMax, tauMinUse, tauMaxUse, taus},
	fMax = Max[freqHz];
	fMin = Min[freqHz];

	tauMin = 1/(2 Pi fMax);
	tauMax = 1/(2 Pi fMin);

	tauMinUse = tauMin * tauMinFactor;
	tauMaxUse = tauMax * tauMaxFactor;

	(* Ensure monotonic valid range *)
	If[tauMinUse <= 0, tauMinUse = tauMin];
	If[tauMaxUse <= tauMinUse, tauMaxUse = tauMax];

	taus = 10.^Range[Log10[tauMinUse], Log10[tauMaxUse], 1./binsPerDecade];
	Developer`ToPackedArray[taus]
];

(* Second-difference matrix for smoothing the distribution (not C0). *)
secondDifferenceMatrix[nBins_Integer] := SparseArray[
	Join[
		Table[{k, k}   ->  1, {k, 1, nBins - 2}],
		Table[{k, k+1} -> -2, {k, 1, nBins - 2}],
		Table[{k, k+2} ->  1, {k, 1, nBins - 2}]
	],
	{nBins - 2, nBins}
];

(* ========================= *)
(* Main API                  *)
(* ========================= *)

Options[DCTSpectrum] = {
	"BinsPerDecade" -> 25,          (* tau grid density *)
	"FMinUse" -> 0.5,
	"FMaxUse" -> 5000,
	"TauMinFactor"  -> 0.1,         (* tauMinUse = tauMin * factor *)
	"TauMaxFactor"  -> 10.0,          (* tauMaxUse = tauMax * factor *)
	"LambdaND"      -> 10^-2,       (* dimensionless regularization strength *)
	"TopPointsForRs"-> 7,           (* how many highest-f points to use for Rs estimate *)
	"MinFreqRsFit" -> minFreqRsFit,
"MaxFreqRsFit" -> maxFreqRsFit
	"WeightMode"    -> "AbsYHalf",   (* choose one of: "AbsYHalf","AbsZHalf","AbsYOne","AbsInvZOne" *)
	"WeightPower"   -> 3/4,          (* wY = |Yint|^(-WeightPower); use 1/2 as a robust default *)
	"Debug" -> False
};
DCTSpectrum[file_String, OptionsPattern[]] := Catch@Module[
  {
    debug,

    (* options *)
    binsPerDecade, tauMinFactor, tauMaxFactor, lambdaND, nTopRs, wPow,
    fMinUse, fMaxUse,

    (* imported data *)
    dat, freqHz, zData, omega, keep, timeS, finishTime,
	
	(*  frequency range for Rs fit *)
	minFreqRsFit, maxFreqRsFit,
    
    (* series resistance + interface admittance *)
    rs, zInt, yInt, wY, wYsqrt,

    (* tau grid + kernels *)
    tauBins, nBins, dLog10, kMat, kUse, aMat, aW, bW, aRI, bRI,

    (* smoothing + NNLS *)
    d2, dFull, aCols, colNormSq, sA2, lambda, dFullN, aAug, bAug, xBest,
    c0, gFit, ciFit,

    (* reconstruct fit *)
    yIntFit, zFit,

    out
  },

  debug = OptionValue["Debug"];
  dbg[debug, "File", file];

  binsPerDecade = OptionValue["BinsPerDecade"];
  tauMinFactor  = OptionValue["TauMinFactor"];
  tauMaxFactor  = OptionValue["TauMaxFactor"];
  lambdaND      = OptionValue["LambdaND"];
  nTopRs        = OptionValue["TopPointsForRs"];
  wPow          = OptionValue["WeightPower"];
  fMinUse       = OptionValue["FMinUse"];
  fMaxUse       = OptionValue["FMaxUse"];
	minFreqRsFit = OptionValue["MinFreqRsFit"];
	maxFreqRsFit = OptionValue["MaxFreqRsFit"];

  (* ---------- Import ---------- *)
  dat    = importEISTxt[file];
  freqHz = dat["FreqHz"];
  zData  = Developer`ToPackedArray[dat["Z"]];
  omega  = Developer`ToPackedArray[2 Pi freqHz];

timeS = dat["TimeS"];
finishTime = Max[timeS];

  dbg[debug, "freqRangeHz BEFORE f-window", {Min[freqHz], Max[freqHz]}];
  With[{n = Min[10, Length[freqHz]]},
    dbg[debug, "highest 10 freqs", Take[Sort[freqHz], -n]];
  ];
  
  
  (* ---------- Rs estimate ---------- *)
  (* rs = estimateRs[freqHz, zData, nTopRs]; *)
  rs = estimateRs[freqHz, zData, minFreqRsFit, maxFreqRsFit];

  dbg[debug, "estimated rs", rs];

  (* ---------- Restrict frequency range ---------- *)
  dbg[debug, "fMinUse", fMinUse];
  dbg[debug, "fMaxUse", fMaxUse];

  keep = (fMinUse <= # <= fMaxUse) & /@ freqHz;
  freqHz = Pick[freqHz, keep];
  zData  = Pick[zData,  keep];
  omega  = Developer`ToPackedArray[2 Pi freqHz];

  dbg[debug, "nPoints after f-window", Length[freqHz]];
  dbg[debug, "freqRangeHz after f-window", {Min[freqHz], Max[freqHz]}];

  assert[Length[freqHz] >= 5, "Too few points after f-window"];
  assert[VectorQ[freqHz, NumericQ], "freqHz not numeric"];
  assert[VectorQ[zData, NumericQ], "zData not numeric"];
  dbg[debug, "firstZ", zData[[1]]];


  zInt = zData - rs;

  (* drop pathological zeros *)
  With[{good = Select[Range[Length[zInt]], Abs[zInt[[#]]] > 0 &]},
    zInt   = zInt[[good]];
    zData  = zData[[good]];
    freqHz = freqHz[[good]];
    omega  = omega[[good]];
  ];

  yInt = Developer`ToPackedArray[1/zInt];

  (* ---------- Weights ---------- *)
  wY = Abs[zData];
  wY = Developer`ToPackedArray[wY];
  wYsqrt = wY^wPow;

  dbg[debug, "WeightPower", wPow];
  dbg[debug, "wYsqrt finite?", FreeQ[wYsqrt, _ComplexInfinity | _DirectedInfinity | Indeterminate | Infinity]];
  assert[VectorQ[wYsqrt, NumericQ], "wYsqrt not numeric"];

  (* ---------- Tau grid ---------- *)
  tauBins = buildTauGrid[freqHz, binsPerDecade, tauMinFactor, tauMaxFactor];
  nBins   = Length[tauBins];
  dLog10  = Developer`ToPackedArray[ConstantArray[1./binsPerDecade, nBins]];

  dbg[debug, "nBins", nBins];
  dbg[debug, "tauRange", {Min[tauBins], Max[tauBins]}];

  (* ---------- Kernel and design matrix ---------- *)
  kMat = Developer`ToPackedArray@Table[
    (I*omega[[j]])/(1 + I*omega[[j]]*tauBins[[k]]),
    {j, Length[omega]}, {k, nBins}
  ];

  kUse = kMat . SparseArray@DiagonalMatrix[dLog10];
  aMat = Join[Transpose[{I*omega}], kUse, 2];

  aW = DiagonalMatrix[wYsqrt] . aMat;
  bW = wYsqrt * yInt;

  aRI = Join[Re[aW], Im[aW]];
  bRI = Join[Re[bW], Im[bW]];

  dbg[debug, "Dimensions(aRI)", Dimensions[aRI]];
  dbg[debug, "Length(bRI)", Length[bRI]];
  assert[MatrixQ[aRI, NumericQ], "aRI not numeric matrix"];
  assert[VectorQ[bRI, NumericQ], "bRI not numeric vector"];
  assert[Dimensions[aRI][[1]] == Length[bRI], "aRI rows != bRI length"];

  (* ---------- Tikhonov smoothing ---------- *)
  d2    = secondDifferenceMatrix[nBins];
  dFull = ArrayFlatten[{{ConstantArray[0., {nBins - 2, 1}], d2}}];

  aCols     = aRI[[All, 2 ;;]];
  colNormSq = Total[aCols^2, {1}];
  sA2       = Median[colNormSq];
  lambda    = lambdaND * sA2;

  dbg[debug, "lambdaND", lambdaND];
  dbg[debug, "sA2", sA2];
  dbg[debug, "lambda", lambda];

  dFullN = N @ Normal @ dFull;

  aAug = Join[aRI, Sqrt[lambda] * dFullN];
  bAug = Join[bRI, ConstantArray[0., nBins - 2]];

  dbg[debug, "Dimensions(aAug)", Dimensions[aAug]];
  dbg[debug, "Length(bAug)", Length[bAug]];
  assert[MatrixQ[aAug, NumericQ], "aAug not numeric matrix"];
  assert[VectorQ[bAug, NumericQ], "bAug not numeric vector"];
  assert[Dimensions[aAug][[1]] == Length[bAug], "aAug rows != bAug length"];

  dbg[debug, "NNLS: A finite?", FreeQ[aAug, _ComplexInfinity | _DirectedInfinity | Indeterminate | Infinity]];
  dbg[debug, "NNLS: b finite?", FreeQ[bAug, _ComplexInfinity | _DirectedInfinity | Indeterminate | Infinity]];

  If[
    !FreeQ[aAug, _ComplexInfinity | _DirectedInfinity | Indeterminate | Infinity] ||
    !FreeQ[bAug, _ComplexInfinity | _DirectedInfinity | Indeterminate | Infinity],
    Throw[Failure["DCTFit", <|"Message" -> "Non-finite values in A or b before NNLS."|>]]
  ];

  (* ---------- NNLS solve (vetted package) ---------- *)
  xBest = Quiet @ Check[NNLSFit`NNLS[aAug, bAug], $Failed];

  dbg[debug, "NNLS: xBest head", Head[xBest]];
  dbg[debug, "NNLS: xBest length", If[ListQ[xBest], Length[xBest], Missing["NotAList"]]];

  If[xBest === $Failed,
    Throw[Failure["DCTFit", <|"Message" -> "NNLSFit`NNLS returned $Failed"|>]]
  ];
  If[!VectorQ[xBest, NumericQ],
    Throw[Failure["DCTFit", <|"Message" -> "NNLSFit`NNLS returned non-numeric xBest"|>]]
  ];
  If[Length[xBest] =!= (nBins + 1),
    Throw[Failure["DCTFit", <|"Message" -> "NNLSFit`NNLS returned wrong-size xBest"|>]]
  ];

  dbg[debug, "NNLS: min(xBest)", Min[xBest]];
  dbg[debug, "NNLS: max(xBest)", Max[xBest]];

  (* ---------- Unpack solution ---------- *)
  c0   = xBest[[1]];
  gFit = xBest[[2 ;;]];
  ciFit = gFit * dLog10;

  (* ---------- Reconstruct fit ---------- *)
  yIntFit = (I*omega)*c0 + kUse . gFit;
  zFit    = rs + 1/yIntFit;

  out =
    <|
      "Tau"      -> tauBins,
      "g"        -> gFit,
      "C0"       -> c0,
      "Rs"       -> rs ,
      "FreqHz"   -> freqHz,
      "ZData"    -> zData,
      "ZFit"     -> zFit,
      "YIntData" -> yInt,
      "YIntFit"  -> yIntFit,
      "FinishTimeS" -> finishTime

    |>;

  dbg[debug, "DCTSpectrum return head", Head[out]];
  dbg[debug, "DCTSpectrum keys", Keys[out]];

  out
];

End[];
EndPackage[];
