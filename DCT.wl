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
Module[{here, nnlsPaths},

  here = DirectoryName[$InputFileName];

  (* search recursively for NNLSFit.m or NNLSFit.wl *)
  nnlsPaths = FileNames[{"NNLSFit.m", "NNLSFit.wl"}, here, Infinity];

  If[nnlsPaths === {},
    Print["NNLSFit not found under: ", here];
    Abort[],
    Get[First[nnlsPaths]]
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


ClearAll[toNumber, splitEISLine, importEISTxt];

toNumber[x_] := Module[{s, y},
  s = StringTrim[ToString[x]];
  If[s === "", Return[Missing["NotNumeric"]]];
  y = Quiet @ Check[ToExpression[s], $Failed];
  If[NumericQ[y], y, Missing["NotNumeric"]]
];

splitEISLine[s_String] := Module[{t = StringTrim[s]},
  Which[
    StringContainsQ[t, ";"],
      StringTrim /@ StringSplit[t, ";"],
    StringContainsQ[t, "\t"],
      StringTrim /@ StringSplit[t, "\t"],
    True,
      DeleteCases[StringSplit[t, Whitespace], ""]
  ]
];

importEISTxt[file_String] := Module[
  {
    rawLines, lines, header, hasHeader,
    splitHeader, colMap,
    freqCol, zreCol, zimCol, timeCol,
    splitRows, numRows, good,
    freqHz, zre, zimNeg, z, timeS, ord
  },

  rawLines = Import[file, "Lines"];
  rawLines = Select[rawLines, StringTrim[#] =!= "" &];

  If[rawLines === {},
    Return[
      Failure["DCTImport", <|"Message" -> "File is empty."|>]
    ]
  ];

  header = First[rawLines];
  hasHeader = StringContainsQ[header, "Frequency", IgnoreCase -> True];

  lines = If[hasHeader, Rest[rawLines], rawLines];

  If[hasHeader,
    splitHeader = splitEISLine[header];
    colMap = AssociationThread[splitHeader -> Range[Length[splitHeader]]];

    freqCol = Lookup[colMap, "Frequency (Hz)", Missing["NoCol"]];
    zreCol  = Lookup[colMap, "Z' (\[CapitalOmega])", Missing["NoCol"]];
    zimCol  = Lookup[colMap, "-Z'' (\[CapitalOmega])", Missing["NoCol"]];
    timeCol = Lookup[colMap, "Time (s)", Missing["NoCol"]];

    If[MemberQ[{freqCol, zreCol, zimCol}, Missing["NoCol"]],
      Return[
        Failure[
          "DCTImport",
          <|
            "Message" -> "Header found but required columns not located.",
            "HeaderFields" -> splitHeader
          |>
        ]
      ]
    ],
    
    (* legacy no-header format *)
    freqCol = 2;
    zreCol  = 3;
    zimCol  = 4;
    timeCol = Missing["NoCol"];
  ];

  splitRows = splitEISLine /@ lines;
  numRows = (toNumber /@ #) & /@ splitRows;

  good = Select[
    numRows,
    Length[#] >= Max[freqCol, zreCol, zimCol] &&
    NumericQ[#[[freqCol]]] &&
    NumericQ[#[[zreCol]]] &&
    NumericQ[#[[zimCol]]] &
  ];

  If[good === {},
    Return[
      Failure["DCTImport", <|"Message" -> "No valid numeric data rows were parsed."|>]
    ]
  ];

  freqHz = good[[All, freqCol]];
  zre    = good[[All, zreCol]];
  zimNeg = good[[All, zimCol]];
  z      = zre + I*(-zimNeg);

  timeS =
    If[
      IntegerQ[timeCol] && Length[First[good]] >= timeCol,
      good[[All, timeCol]],
      ConstantArray[Missing["NoTime"], Length[good]]
    ];

  ord = Ordering[freqHz];

  <|
    "FreqHz" -> freqHz[[ord]],
    "Z" -> z[[ord]],
    "TimeS" -> timeS[[ord]]
  |>
];


(* Estimate Rs from points in a specified frequency window.
   Used only as an INITIAL GUESS for the OUTER LOOP. *)
estimateRs[freqHz_List, z_List, minFreqRsFit_?NumericQ, maxFreqRsFit_?NumericQ] := Module[
	{keep, top, rs},
	keep = (minFreqRsFit <= # <= maxFreqRsFit) & /@ freqHz;
	top = Re[Pick[z, keep]];
	top = Select[top, NumericQ];
	If[top === {}, Return[0.]];
	(* rs = Median[top]; *)
	rs = Min[top];
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

	If[tauMinUse <= 0, tauMinUse = tauMin];
	If[tauMaxUse <= tauMinUse, tauMaxUse = tauMax];

	taus = 10.^Range[Log10[tauMinUse], Log10[tauMaxUse], 1./binsPerDecade];
	Developer`ToPackedArray[taus]
];

(* Second-difference matrix for smoothing the distribution (not C0). *)
secondDifferenceMatrix[nBins_Integer] := SparseArray[
	Join[
		Table[{k, k} -> 1, {k, 1, nBins - 2}],
		Table[{k, k + 1} -> -2, {k, 1, nBins - 2}],
		Table[{k, k + 2} -> 1, {k, 1, nBins - 2}]
	],
	{nBins - 2, nBins}
];

(* ========================= *)
(* Main API                  *)
(* ========================= *)

Options[DCTSpectrum] = {
	"BinsPerDecade" -> 25,
	"FMinUse" -> 0.5,
	"FMaxUse" -> 5000,
	"TauMinFactor" -> 0.1,
	"TauMaxFactor" -> 10.0,
	"LambdaND" -> 10^-2,
	"TopPointsForRs" -> 7,
	"WeightMode" -> "AbsYHalf",
	"WeightPower" -> 3/4,
	"Debug" -> False
};

DCTSpectrum[file_String, OptionsPattern[]] := Catch@Module[
	{
		debug,

		(* options *)
		binsPerDecade, tauMinFactor, tauMaxFactor, lambdaND, nTopRs, wPow,
		fMinUse, fMaxUse,

		(* imported data *)
		dat, freqHzAll, zDataAll, timeS, finishTime,

		(* INNER LOOP data: DCT fit range *)
		keepDCT, freqHzDCT, zDataDCT, omegaDCT,

		(* OUTER LOOP data: Rs-fit range *)
		keepRs, freqHzRs, zDataRs, omegaRs,

		(* OUTER LOOP quantities *)
		rs0, rsLower, rsUpper, rsBest, rsWindowRe,
		outerObjective, bestSolve,

		(* INNER LOOP helper *)
		solveLadderGivenRs,

		out
	},

	debug = OptionValue["Debug"];
	dbg[debug, "File", file];

	binsPerDecade = OptionValue["BinsPerDecade"];
	tauMinFactor = OptionValue["TauMinFactor"];
	tauMaxFactor = OptionValue["TauMaxFactor"];
	lambdaND = OptionValue["LambdaND"];
	nTopRs = OptionValue["TopPointsForRs"];
	wPow = OptionValue["WeightPower"];
	fMinUse = N[OptionValue["FMinUse"]];
	fMaxUse = N[OptionValue["FMaxUse"]];

	(* ---------- Import full data ---------- *)
	dat = importEISTxt[file];
	freqHzAll = dat["FreqHz"];
	zDataAll = Developer`ToPackedArray[dat["Z"]];

	timeS = dat["TimeS"];
	finishTime = Max[timeS];

	dbg[debug, "full freq range", {Min[freqHzAll], Max[freqHzAll]}];

	(* ======================================================== *)
	(* Build INNER LOOP data: DCT fit range                     *)
	(* Maxwell ladder admittance is fit ONLY on this range      *)
	(* ======================================================== *)
	keepDCT = (fMinUse <= # <= fMaxUse) & /@ freqHzAll;
	freqHzDCT = Pick[freqHzAll, keepDCT];
	zDataDCT = Pick[zDataAll, keepDCT];
	omegaDCT = Developer`ToPackedArray[2 Pi freqHzDCT];

	dbg[debug, "INNER LOOP fMinUse", fMinUse];
	dbg[debug, "INNER LOOP fMaxUse", fMaxUse];
	dbg[debug, "INNER LOOP nPoints", Length[freqHzDCT]];
	dbg[debug, "INNER LOOP freq range", {Min[freqHzDCT], Max[freqHzDCT]}];

	assert[Length[freqHzDCT] >= 5, "Too few points in DCT fit range"];
	assert[VectorQ[freqHzDCT, NumericQ], "freqHzDCT not numeric"];
	assert[VectorQ[zDataDCT, NumericQ], "zDataDCT not numeric"];

	(* ======================================================== *)
	(* Build OUTER LOOP data: Rs fit range                      *)
	(* impedance residual for Rs is evaluated ONLY on this      *)
	(* separate frequency window                                *)
	(* ======================================================== *)
	keepRs = (fMinUse <= # <= fMaxUse) & /@ freqHzAll;
	freqHzRs = Pick[freqHzAll, keepRs];
	zDataRs = Pick[zDataAll, keepRs];
	omegaRs = Developer`ToPackedArray[2 Pi freqHzRs];

	dbg[debug, "OUTER LOOP nPoints", Length[freqHzRs]];
	dbg[debug, "OUTER LOOP freq range", If[Length[freqHzRs] > 0, {Min[freqHzRs], Max[freqHzRs]}, Missing["NoPoints"]]];

	assert[Length[freqHzRs] >= 1, "Too few points in Rs fit range"];
	assert[VectorQ[freqHzRs, NumericQ], "freqHzRs not numeric"];
	assert[VectorQ[zDataRs, NumericQ], "zDataRs not numeric"];

	(* ======================================================== *)
	(* INNER LOOP: solve Maxwell ladder in admittance space     *)
	(* for a FIXED candidate Rs                                 *)
	(* Uses ONLY DCT range and ONLY DCT tau grid                *)
	(* ======================================================== *)
	solveLadderGivenRs[rsCand_?NumericQ] := Module[
		{
			zIntDCT, yIntDCT, yIntDataRs, wY, wYsqrt, wYRs,
			tauBins, nBins, dLog10,
			kMatDCT, kUseDCT, aMat, aW, bW, aRI, bRI,
			d2, dFull, aCols, colNormSq, sA2, lambda, dFullN, aAug, bAug, xBest,
			c0, gFit, ciFit, yIntFitDCT, zFitDCT,
			kMatRs, kUseRs, yIntFitRs,
			objZ, objY
		},

		zIntDCT = zDataDCT - rsCand;

		If[Min[Abs[zIntDCT]] <= 10^-15,
			Return[<|"OK" -> False, "Message" -> "Z - Rs too small in DCT range", "ObjZ" -> Infinity|>]
		];

(*
		yIntDCT = Developer`ToPackedArray[1/zIntDCT];

		(* ---------- Weights on DCT range ---------- *)
		(* wY = Developer`ToPackedArray[Abs[zDataDCT]]; *)
		wY = Developer`ToPackedArray[Abs[zIntDCT]];
		wYsqrt = wY^wPow;
*)		

		yIntDCT = Developer`ToPackedArray[1/zIntDCT];
		wYsqrt = Developer`ToPackedArray[Abs[yIntDCT]^wPow];

		If[!VectorQ[wYsqrt, NumericQ],
			Return[<|"OK" -> False, "Message" -> "wYsqrt not numeric", "ObjZ" -> Infinity|>]
		];
		If[!FreeQ[wYsqrt, _ComplexInfinity | _DirectedInfinity | Indeterminate | Infinity],
			Return[<|"OK" -> False, "Message" -> "wYsqrt non-finite", "ObjZ" -> Infinity|>]
		];

		(* ---------- Tau grid built ONLY from DCT range ---------- *)
		tauBins = buildTauGrid[freqHzDCT, binsPerDecade, tauMinFactor, tauMaxFactor];
		nBins = Length[tauBins];
		dLog10 = Developer`ToPackedArray[ConstantArray[1./binsPerDecade, nBins]];

		(* ---------- Kernel and design matrix on DCT range ---------- *)
		kMatDCT = Developer`ToPackedArray@Table[
			(I*omegaDCT[[j]])/(1 + I*omegaDCT[[j]]*tauBins[[k]]),
			{j, Length[omegaDCT]}, {k, nBins}
		];

		kUseDCT = kMatDCT . SparseArray@DiagonalMatrix[dLog10];
		aMat = Join[Transpose[{I*omegaDCT}], kUseDCT, 2];

		aW = DiagonalMatrix[wYsqrt] . aMat;
		bW = wYsqrt * yIntDCT;

		aRI = Join[Re[aW], Im[aW]];
		bRI = Join[Re[bW], Im[bW]];

		If[!MatrixQ[aRI, NumericQ] || !VectorQ[bRI, NumericQ],
			Return[<|"OK" -> False, "Message" -> "aRI/bRI not numeric", "ObjZ" -> Infinity|>]
		];
		If[Dimensions[aRI][[1]] =!= Length[bRI],
			Return[<|"OK" -> False, "Message" -> "aRI rows != bRI length", "ObjZ" -> Infinity|>]
		];

		(* ---------- Tikhonov smoothing on ladder only ---------- *)
		d2 = secondDifferenceMatrix[nBins];
		dFull = ArrayFlatten[{{ConstantArray[0., {nBins - 2, 1}], d2}}];

		aCols = aRI[[All, 2 ;;]];
		colNormSq = Total[aCols^2, {1}];
		sA2 = Median[colNormSq];
		lambda = lambdaND * sA2;

		dFullN = N @ Normal @ dFull;
		aAug = Join[aRI, Sqrt[lambda] * dFullN, 1];
		bAug = Join[bRI, ConstantArray[0., nBins - 2]];

		If[
			!MatrixQ[aAug, NumericQ] || !VectorQ[bAug, NumericQ] ||
			Dimensions[aAug][[1]] =!= Length[bAug],
			Return[<|"OK" -> False, "Message" -> "aAug/bAug invalid", "ObjZ" -> Infinity|>]
		];

		If[
			!FreeQ[aAug, _ComplexInfinity | _DirectedInfinity | Indeterminate | Infinity] ||
			!FreeQ[bAug, _ComplexInfinity | _DirectedInfinity | Indeterminate | Infinity],
			Return[<|"OK" -> False, "Message" -> "aAug/bAug non-finite", "ObjZ" -> Infinity|>]
		];

		(* ---------- NNLS solve ---------- *)
		xBest = Quiet @ Check[NNLSFit`NNLS[aAug, bAug], $Failed];

		If[xBest === $Failed,
			Return[<|"OK" -> False, "Message" -> "NNLS returned $Failed", "ObjZ" -> Infinity|>]
		];
		If[!VectorQ[xBest, NumericQ],
			Return[<|"OK" -> False, "Message" -> "NNLS returned non-numeric", "ObjZ" -> Infinity|>]
		];
		If[Length[xBest] =!= (nBins + 1),
			Return[<|"OK" -> False, "Message" -> "NNLS returned wrong-size xBest", "ObjZ" -> Infinity|>]
		];

		c0 = xBest[[1]];
		gFit = xBest[[2 ;;]];
		ciFit = gFit * dLog10;

		(* ---------- Reconstruct fit on DCT range ---------- *)
		yIntFitDCT = (I*omegaDCT)*c0 + kUseDCT . gFit;
		If[Min[Abs[yIntFitDCT]] <= 10^-30,
			Return[<|"OK" -> False, "Message" -> "YintFit too small", "ObjZ" -> Infinity|>]
		];

		zFitDCT = rsCand + 1/yIntFitDCT;

		(* ---------- Reconstruct SAME ladder on Rs-fit range ---------- *)
		kMatRs = Developer`ToPackedArray@Table[
			(I*omegaRs[[j]])/(1 + I*omegaRs[[j]]*tauBins[[k]]),
			{j, Length[omegaRs]}, {k, nBins}
		];

		kUseRs = kMatRs . SparseArray@DiagonalMatrix[dLog10];
		yIntFitRs = (I*omegaRs)*c0 + kUseRs . gFit;

		If[Min[Abs[yIntFitRs]] <= 10^-30,
			Return[<|"OK" -> False, "Message" -> "YintFitRs too small", "ObjZ" -> Infinity|>]
		];

		yIntDataRs = Developer`ToPackedArray[1/(zDataRs - rsCand)];

		If[!VectorQ[yIntDataRs, NumericQ] ||
		   !FreeQ[yIntDataRs, _ComplexInfinity | _DirectedInfinity | Indeterminate | Infinity],
		   Return[<|"OK" -> False, "Message" -> "bad outer-loop admittance data", "ObjZ" -> Infinity|>]
		];
		
		wYRs = Developer`ToPackedArray[Abs[yIntDataRs]^wPow];
		
		objZ =
			Total[
				wYRs *
				(
					Re[yIntDataRs - yIntFitRs]^2 +
					Im[yIntDataRs - yIntFitRs]^2
				)
			];
			   

		(* ---------- diagnostic objective in admittance space, ONLY on DCT range ---------- *)
		objY = Total[Re[yIntDCT - yIntFitDCT]^2 + Im[yIntDCT - yIntFitDCT]^2];

		<|
			"OK" -> True,
			"Message" -> "",
			"Rs" -> rsCand,
			"Tau" -> tauBins,
			"g" -> gFit,
			"C0" -> c0,
			"FreqHz" -> freqHzDCT,
			"ZData" -> zDataDCT,
			"ZFit" -> zFitDCT,
			"YIntData" -> yIntDCT,
			"YIntFit" -> yIntFitDCT,
			"FinishTimeS" -> finishTime,
			"ObjZ" -> objZ,
			"ObjY" -> objY
		|>
	];

	(* ======================================================== *)
	(* OUTER LOOP: solve for Rs using ONLY Rs-fit frequency band *)
	(* ======================================================== *)

	rs0 = estimateRs[freqHzAll, zDataAll, fMinUse, fMaxUse]; 
	If[!NumericQ[rs0], rs0 = Max[0., Min[Re[zDataRs]]]];
	dbg[debug, "initial Rs guess", rs0];

	rsWindowRe = Re[zDataRs];
	rsWindowRe = Select[rsWindowRe, NumericQ];

	rsLower = Max[0., Min[rsWindowRe] - 0.5 (Max[rsWindowRe] - Min[rsWindowRe])];
	rsUpper = Max[rsLower + 1., Max[rsWindowRe] + 0.5 (Max[rsWindowRe] - Min[rsWindowRe])];

	dbg[debug, "OUTER LOOP Rs bounds", {rsLower, rsUpper}];

	Module[
		{
			callCount = 0,
			rsCache = <||>,
			rsKey,
			t0, sol, objVal,
			refineSol
		},

		outerObjective[rsVar_?NumericQ] := Module[{},

			callCount++;
			rsKey = ToString @ NumberForm[N[rsVar], {16, 6}];

			If[KeyExistsQ[rsCache, rsKey],
				If[debug,
					Print[
						"[OUTER DEBUG] cached call ", callCount,
						"   Rs=", N[rsVar],
						"   Obj=", rsCache[rsKey]
					];
				];
				Return[rsCache[rsKey]];
			];

			t0 = AbsoluteTime[];
			sol = solveLadderGivenRs[rsVar];

			objVal =
				If[TrueQ[sol["OK"]],
					N[sol["ObjZ"]],
					Infinity
				];

			rsCache[rsKey] = objVal;

			If[debug,
				Print[
					"[OUTER DEBUG] call ", callCount,
					"   Rs=", N[rsVar],
					"   Obj=", objVal,
					"   dt=", NumberForm[AbsoluteTime[] - t0, {6, 3}], " s",
					"   cache size=", Length[rsCache]
				];
			];

			objVal
		];

		dbg[debug, "OUTER LOOP objective at initial Rs", outerObjective[rs0]];

		refineSol = Quiet @ Check[
			FindMinimum[
				{
					outerObjective[r],
					rsLower <= r <= rsUpper
				},
				{{r, rs0}},
				Method -> "Brent",
				MaxIterations -> 40,
				WorkingPrecision -> MachinePrecision
			],
			$Failed
		];

		If[refineSol === $Failed || !MatchQ[refineSol, {_?NumericQ, {__Rule}}],
			rsBest = rs0;,
			rsBest = r /. refineSol[[2]]
		];

		If[!NumericQ[rsBest], rsBest = rs0];
		rsBest = N[rsBest];

		dbg[debug, "OUTER LOOP best Rs", rsBest];
		dbg[debug, "OUTER LOOP best objective", outerObjective[rsBest]];
		dbg[debug, "OUTER LOOP total objective calls", callCount];
		dbg[debug, "OUTER LOOP unique cached Rs values", Length[rsCache]];
	];

	If[rsBest === $Failed || !NumericQ[rsBest],
		Throw[Failure["DCTFit", <|"Message" -> "OUTER LOOP Rs optimization failed."|>]]
	];

	bestSolve = solveLadderGivenRs[rsBest];

	If[!TrueQ[bestSolve["OK"]],
		Throw[Failure["DCTFit", <|"Message" -> "INNER LOOP solve failed at optimized Rs."|>]]
	];

	out =
		<|
			"Tau" -> bestSolve["Tau"],
			"g" -> bestSolve["g"],
			"C0" -> bestSolve["C0"],
			"Rs" -> bestSolve["Rs"],
			"FreqHz" -> bestSolve["FreqHz"],
			"ZData" -> bestSolve["ZData"],
			"ZFit" -> bestSolve["ZFit"],
			"YIntData" -> bestSolve["YIntData"],
			"YIntFit" -> bestSolve["YIntFit"],
			"FinishTimeS" -> bestSolve["FinishTimeS"],
			"ObjZ" -> bestSolve["ObjZ"],
			"ObjY" -> bestSolve["ObjY"],
			"RsFitFreqHz" -> freqHzRs,
			"RsFitZData" -> zDataRs
		|>;

	dbg[debug, "DCTSpectrum return head", Head[out]];
	dbg[debug, "DCTSpectrum keys", Keys[out]];

	out
];
End[];
EndPackage[];
