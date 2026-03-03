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
	freqCol, zreCol, zimCol, splitRows, numRows, good, freqHz, zre, zimNeg, z},

	rawLines = Import[file, "Lines"];
	rawLines = Select[rawLines, StringTrim[#] =!= "" &];

	(* Detect header: first line contains "Frequency" *)
	header = First[rawLines];
	hasHeader = StringContainsQ[header, "Frequency"];

	lines = If[hasHeader, Rest[rawLines], rawLines];

	(* If we have a header, locate columns by name; else fall back to old positions *)
	If[hasHeader,
		splitHeader = StringTrim /@ StringSplit[header, ";"];
		colMap = AssociationThread[splitHeader -> Range[Length[splitHeader]]];

		(* allow a couple common variants *)
		freqCol = Lookup[colMap, "Frequency (Hz)", Missing["NoCol"]];
		zreCol  = Lookup[colMap, "Z' (\[CapitalOmega])", Missing["NoCol"]];
		zimCol  = Lookup[colMap, "-Z'' (\[CapitalOmega])", Missing["NoCol"]];

		If[MemberQ[{freqCol, zreCol, zimCol}, Missing["NoCol"]],
			Throw[Failure["DCTImport", <|"Message" -> "Header found but required columns not located."|>]]
		],
		(* no header: assume legacy positions *)
		freqCol = 2; zreCol = 3; zimCol = 4;
	];

	splitRows = StringSplit[#, ";"] & /@ lines;
	numRows = (toNumber /@ #) & /@ splitRows;

	good = Select[numRows, (Length[#] >= Max[freqCol, zreCol, zimCol] && 
		NumericQ[#[[freqCol]]] && NumericQ[#[[zreCol]]] && NumericQ[#[[zimCol]]]) &];

	freqHz = good[[All, freqCol]];
	zre    = good[[All, zreCol]];
	zimNeg = good[[All, zimCol]];

	z = zre + I*(-zimNeg);  (* file stores -Zim *)

	With[{ord = Ordering[freqHz]},
		<|"FreqHz" -> freqHz[[ord]], "Z" -> z[[ord]]|>
	]
];

(* Estimate Rs from the highest-frequency points.
   This is intentionally simple: take the median of Re[Z] over the top N points. *)
estimateRs[freqHz_List, z_List, nTop_Integer] := Module[
	{n = Length[freqHz], ord, top, rs},
	ord = Ordering[freqHz, -Min[nTop, n]];
	top = Re[z[[ord]]];
	rs = Median[top];
	If[NumericQ[rs] && rs >= 0, rs, 0.]
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
	"TauMinFactor"  -> 10.,         (* tauMinUse = tauMin * factor *)
	"TauMaxFactor"  -> 1.,          (* tauMaxUse = tauMax * factor *)
	"LambdaND"      -> 10^-2,       (* dimensionless regularization strength *)
	"TopPointsForRs"-> 7,           (* how many highest-f points to use for Rs estimate *)
	"WeightMode"    -> "AbsYHalf",   (* choose one of: "AbsYHalf","AbsZHalf","AbsYOne","AbsInvZOne" *)
	"WeightPower"   -> 3/4,          (* wY = |Yint|^(-WeightPower); use 1/2 as a robust default *)
	"Debug" -> False
};

DCTSpectrum[file_String, OptionsPattern[]] := Catch@Module[
  {
    debug,

    (* options *)
    binsPerDecade, tauMinFactor, tauMaxFactor, lambdaND, nTopRs, wPow,

    (* imported data *)
    dat, freqHz, zData, omega,

    (* series resistance + interface admittance *)
    rs, zInt, yInt, wY, wYsqrt,

    (* tau grid + kernels *)
    tauBins, nBins, dLog10, kMat, kUse, aMat, aW, bW, aRI, bRI,

    (* smoothing + NNLS *)
    d2, dFull, aAug, bAug, xVars, obj, sol, xBest, c0, gFit, ciFit,

    (* reconstruct fit *)
    yIntFit, zFit
  },

	debug = OptionValue["Debug"];
	dbg[debug, "File", file];
	
	binsPerDecade = OptionValue["BinsPerDecade"];
	tauMinFactor  = OptionValue["TauMinFactor"];
	tauMaxFactor  = OptionValue["TauMaxFactor"];
	lambdaND      = OptionValue["LambdaND"];
	nTopRs        = OptionValue["TopPointsForRs"];
	wPow          = OptionValue["WeightPower"];
	
	dat    = importEISTxt[file];
	freqHz = dat["FreqHz"];
	zData  = Developer`ToPackedArray[dat["Z"]];
	omega  = Developer`ToPackedArray[2 Pi freqHz];
	dbg[debug, "freqRangeHz BEFORE f-window", {Min[freqHz], Max[freqHz]}];
	
	With[{n = Min[10, Length[freqHz]]},
	 dbg[debug, "highest 10 freqs", Take[Sort[freqHz], -n]];
	];
	  
	(*----restrict frequency range----*)
	fMinUse = OptionValue["FMinUse"];
	fMaxUse = OptionValue["FMaxUse"];
	dbg[debug, "fMinUse", fMinUse];
	  dbg[debug, "fMaxUse", fMaxUse];
	keep = (fMinUse <= # <= fMaxUse) & /@ freqHz;
	
	freqHz = Pick[freqHz, keep];
	zData  = Pick[zData,  keep];
	
	omega  = Developer`ToPackedArray[2 Pi freqHz];
	
	dbg[debug, "nPoints after f-window", Length[freqHz]];
	dbg[debug, "freqRangeHz after f-window", {Min[freqHz], Max[freqHz]}];
	  
	  
	  
	dbg[debug, "nPoints", Length[freqHz]];
	dbg[debug, "freqRangeHz", {Min[freqHz], Max[freqHz]}];
	dbg[debug, "firstZ", zData[[1]]];
	assert[VectorQ[freqHz, NumericQ], "freqHz not numeric"];
	assert[VectorQ[zData, NumericQ], "zData not numeric"];

	rs   = estimateRs[freqHz, zData, nTopRs];
  
	dbg[debug, "estimated rs ", rs];
  
	zInt = zData - rs;

	(* Avoid division by zero: drop any pathological points (rare, but keeps the core solver robust). *)
	With[{good = Select[Range[Length[zInt]], Abs[zInt[[#]]] > 0 &]},
		zInt   = zInt[[good]];
		zData  = zData[[good]];
		freqHz = freqHz[[good]];
		omega  = omega[[good]];
	];

	yInt = Developer`ToPackedArray[1/zInt];

	wY = Abs[zData];  (* matches: wY = Abs[1/Zdata]^(-1) *)
	wY = Developer`ToPackedArray[wY];
	wYsqrt = Sqrt[wY];
	wYsqrt = wY^wPow; (* don't use sqrt weighting!! *)
	
	tauBins = buildTauGrid[freqHz, binsPerDecade, tauMinFactor, tauMaxFactor];
	nBins   = Length[tauBins];
	dLog10  = Developer`ToPackedArray[ConstantArray[1./binsPerDecade, nBins]];

	(* Kernel: (i \[Omega])/(1 + i \[Omega] \[Tau])  for a Maxwell/DCT ladder written in admittance form *)
	kMat = Developer`ToPackedArray@Table[
		(I*omega[[j]])/(1 + I*omega[[j]]*tauBins[[k]]),
		{j, Length[omega]}, {k, nBins}
	];

	(* Convention: solve for g(\[Tau]) in units of F/decade.
     Implemented by scaling columns by \[CapitalDelta]log10 and solving for g directly. *)
	kUse = kMat . SparseArray@DiagonalMatrix[dLog10];

	(* Design matrix includes C0 (the i \[Omega] C0 term) plus the distribution columns. *)
	aMat = Join[Transpose[{I*omega}], kUse, 2];  (* N x (1+nBins) *)

	(* Weight in complex space, then split into real/imag stacked system. *)
	aW = DiagonalMatrix[wYsqrt] . aMat;
	bW = wYsqrt * yInt;

	aRI = Join[Re[aW], Im[aW]];
	bRI = Join[Re[bW], Im[bW]];
	dbg[debug, "Dimensions(aRI)", Dimensions[aRI]];
	dbg[debug, "Length(bRI)", Length[bRI]];
	dbg[debug, "MatrixQ(aRI)", MatrixQ[aRI, NumericQ]];
	dbg[debug, "VectorQ(bRI)", VectorQ[bRI, NumericQ]];
	assert[MatrixQ[aRI, NumericQ], "aRI is not a numeric matrix"];
	assert[VectorQ[bRI, NumericQ], "bRI is not a numeric vector"];
	assert[Dimensions[aRI][[1]] == Length[bRI], "aRI rows != bRI length"];
	

    (* Tikhonov smoothing on the distribution only (not on C0). *)
	d2 = secondDifferenceMatrix[nBins];
	dFull = ArrayFlatten[{{ConstantArray[0., {nBins - 2, 1}], d2}}];

	(* Auto-scale lambda so it is comparable to the data misfit scale. *)
	Module[{aCols, colNormSq, sA2, lambda, dFullN, aAugN, bAugN, Qmat, cvec},

    aCols = aRI[[All, 2 ;;]];                 (* exclude C0 column *)
    colNormSq = Total[aCols^2, {1}];
    sA2 = Median[colNormSq];
    lambda = lambdaND * sA2;

    (* Critical: convert SparseArray to a normal numeric matrix BEFORE Join *)
    dFullN = N @ Normal @ dFull;

    aAugN = Join[aRI, Sqrt[lambda] * dFullN];
    bAugN = Join[bRI, ConstantArray[0., nBins - 2]];

	dbg[debug, "Dimensions(aAugN)", Dimensions[aAugN]];
	dbg[debug, "Length(bAugN)", Length[bAugN]];
	assert[MatrixQ[aAugN, NumericQ], "aAugN not numeric matrix"];
	assert[VectorQ[bAugN, NumericQ], "bAugN not numeric vector"];
	assert[Dimensions[aAugN][[1]] == Length[bAugN], "aAugN rows != bAugN length"];

    (* NNLS via convex QP: minimize ||A x - b||^2 subject to x >= 0 *)
    xVars = Array[x, nBins + 1];

    Qmat = 2.0 * Transpose[aAugN] . aAugN;
    cvec = -2.0 * Transpose[aAugN] . bAugN;

    (* Force exact shapes: Qmat must be 2D, cvec must be 1D *)
    Qmat = N @ Qmat;
    cvec = N @ Flatten[cvec];
	dbg[debug, "Head(Qmat)", Head[Qmat]];
	dbg[debug, "Head(cvec)", Head[cvec]];
	dbg[debug, "Dimensions(Qmat)", Dimensions[Qmat]];
	dbg[debug, "Dimensions(cvec)", Dimensions[cvec]];
	dbg[debug, "Length(xVars)", Length[xVars]];
	dbg[debug, "Qmat[[1,1]] NumericQ", Quiet@NumericQ[Qmat[[1, 1]]]];
	dbg[debug, "cvec[[1]] NumericQ", Quiet@NumericQ[cvec[[1]]]];
	
	(* This catches the exact triple-brace mistake *)
	dbg[debug, "ObjectiveCandidateHeads", {Head[{Qmat, cvec}], Head[Qmat], Head[cvec]}];
	
	assert[MatrixQ[Qmat, NumericQ], "Qmat is not a numeric matrix"];
	assert[VectorQ[cvec, NumericQ], "cvec is not a numeric vector"];
	assert[Dimensions[Qmat] == {Length[xVars], Length[xVars]}, "Qmat shape mismatch"];
	assert[Length[cvec] == Length[xVars], "cvec length mismatch"];

	    (* NNLS via NMinimize *)
    xVars = Array[x, nBins + 1];
	
	dbg[debug, "aAugN finite?", FreeQ[aAugN, _ComplexInfinity | _DirectedInfinity | Indeterminate | Infinity]];
	dbg[debug, "bAugN finite?", FreeQ[bAugN, _ComplexInfinity | _DirectedInfinity | Indeterminate | Infinity]];
	
	If[
	  !FreeQ[aAugN, _ComplexInfinity | _DirectedInfinity | Indeterminate | Infinity] ||
	  !FreeQ[bAugN, _ComplexInfinity | _DirectedInfinity | Indeterminate | Infinity],
	  Throw[Failure["DCTFit", <|"Message" -> "Non-finite values detected in A or b (Inf/Indeterminate). Check weights and Yint."|>]]
	];

    obj = (aAugN . xVars - bAugN) . (aAugN . xVars - bAugN);
	(* ---------- Positivity constraints ---------- *)
	(* Enforce g >= 0 (distribution only). If you also want C0 >= 0, keep the first line. *)

	constraints = Join[
	  {xVars[[1]] >= 0},                 (* C0 >= 0  , remove this line if you want C0 free *)
	  Thread[xVars[[2 ;;]] >= 0]         (* g_k >= 0 for all bins *)
	];


	sol = NMinimize[
	  {obj, constraints},
	  xVars,
	  Method -> Automatic,
	  MaxIterations -> 2000
	];
	
	dbg[debug, "NMinimize result head", Head[sol]];
	dbg[debug, "NMinimize result", sol];

    xBest = xVars /. sol[[2]];
    ];
  
	c0   = xBest[[1]];
	gFit = xBest[[2 ;;]];               (* already F/decade by construction *)
	ciFit = gFit * dLog10;              (* per-bin capacitance, if you ever want it *)

	(* Reconstruct fit in admittance and impedance space. *)
	yIntFit = (I*omega)*c0 + kUse . gFit;
	zFit    = rs + 1/yIntFit;

	<|
		"Tau"    -> tauBins,
		"g"      -> gFit,
		"C0"     -> c0,
		"Rs"     -> rs,
		"FreqHz" -> freqHz,
		"ZData"  -> zData,
		"ZFit"   -> zFit,
		"YIntData" -> yInt,
		"YIntFit"  -> yIntFit    
	|>
];

End[];
EndPackage[];
