(* ::Package:: *)

(* ============================================================ *)
(* DCT batch runner                                              *)
(*   1) Find all .txt files in a directory                        *)
(*   2) Fit each file to a DCT spectrum                            *)
(*   3) Plot the results                                           *)
(* ============================================================ *)


(* ---------- LOAD PACKAGE ---------- *)
ClearAll["DCT`*"];  (* unload package *)

dctPaths = FileNames["DCT.wl", NotebookDirectory[], Infinity];

If[dctPaths === {},
  Print["DCT.wl not found under: ", NotebookDirectory[]];
  Abort[],
  Get[First[dctPaths]]
];


ClearAll["NNLS`*"];  (* unload *)

NNLSPackagePath =
  First @ FileNames["NNLSFit.m", NotebookDirectory[], Infinity];

Get[NNLSPackagePath]


(* ---------- USER SETTINGS ---------- *)












(* ::InheritFromParent:: *)
(**)





(* dataDir = "C:\\Users\\bobla\\Documents\\DCT\\data";
dataDir = "C:\\Users\\bobla\\Documents\\DCT\\data\\2026-02-25-titration";
dataDir = "C:\\Users\\bobla\\Documents\\DCT\\data\\2026-02-25";
dataDir = "C:\\Users\\bobla\\Documents\\DCT\\data\\drift_03252026\\E1"
*)
baseDir = FileNameJoin[{NotebookDirectory[], "data"}] // ExpandFileName


(* examples of switching datasets *)
dataDir = FileNameJoin[{baseDir, "2026-02-25-titration"}];
dataDir = FileNameJoin[{baseDir, "2026-02-25"}];
dataDir = FileNameJoin[{baseDir, "drift_03252026", "E1"}]

debugFlag = False;
fileDecimation = 1;   (* keep every Nth file: 10 -> ~450/10 = 45 files *)

lambdaND = 1 10^-3;

constantPhaseElementFlag = False;

binsPerDecade = 35;


weightPower = -0.5; (* weight each point inverse to variance *)

(* weightPowerResistance = -3.5; *)
weightPowerResistance = weightPower + 4;

fMinUse = 0.5;      (* Hz *) (* starting to see resistive behavior at low freq (oxygen reduction? diffusion?) *)
fMaxUse = 200;     (* Hz *)

minFreqRsFit = 200;  (* Hz *)
maxFreqRsFit = 350;  (* Hz *)



(* Optional: restrict the frequency range (must match package options) *)
paddingDecades = 0.0; (* sets how many decades beyond the measured frequency range the tau values extend *)
(* 0 = no padding, 1 = a decdade of apadding. *)

tauMinFactor = 10^-paddingDecades;
tauMaxFactor = 10^paddingDecades;          (* tauMaxUse = tauMax * factor *) (* add factor of 10 to the maximum as padding *)
(* 
minFreqRsFit = fMaxUse 10^paddingDecades;   (* Hz *)
maxFreqRsFit = 300 + minFreqRsFit;  (* Hz *)
*)


(* ============================================================ *)
(* 1) FIND ALL TXT FILES                                        *)
(* ============================================================ *)

txtFilesAll =
  SortBy[
    FileNames["*.txt", dataDir],
    ToExpression @ First @ StringCases[FileNameTake[#], "(" ~~ x : NumberString ~~ ")" :> x] &
  ];

Print["Found ", Length[txtFilesAll], " .txt files total."];
If[Length[txtFilesAll] == 0, Abort[]];

fileDecimation = Max[1, Round[fileDecimation]];

txtFiles =
  With[{idx = Range[1, Length[txtFilesAll], fileDecimation]},
    txtFilesAll[[idx]]
  ];

Print[
  "Using ", Length[txtFiles],
  " files after decimation (every ", fileDecimation, "th file)."
];

Print["Files used:"];
Print /@ txtFiles;


(* ============================================================ *)
(* 2) FIT EACH FILE TO A DCT SPECTRUM                            *)
(*    - sort by numeric index in parentheses: ... (n).txt        *)
(* ============================================================ *)

(* ---------- helper: "natural" sort by trailing (...) integer ---------- *)
(* Example: "E1_EIS_15(1).txt" < "E1_EIS_15(26).txt" *)
fileOrderKey[path_String] := Module[{base, n},
  base = FileBaseName[path];
  n = Quiet @ Check[
    ToExpression @ StringReplace[base, RegularExpression[".*\\((\\d+)\\)$"] -> "$1"],
    Missing["NoIndex"]
  ];
  If[IntegerQ[n], {StringReplace[base, RegularExpression["\\(\\d+\\)$"] -> ""], n}, {base, Infinity}]
];

(* Sort the decimated list by the numeric index if present *)
txtFiles = SortBy[txtFiles, fileOrderKey];
(* ============================================================ *)
(* 2) FIT EACH FILE TO A DCT SPECTRUM  (with progress indicator) *)
(* ============================================================ *)

nFiles = Length[txtFiles];
i = 0;
currentFile = "";

results = Monitor[
  Table[
    Module[{file = f, spec, msg = "", ok = True},

      i++;
      currentFile = FileNameTake[file];

      spec = Quiet[
        Check[
          DCT`DCTSpectrum[
            file,
            "Debug" -> debugFlag,
            "LambdaND" -> lambdaND,
            "FMinUse" -> fMinUse,
            "FMaxUse" -> fMaxUse,
            "TauMinFactor" -> tauMinFactor,
            "TauMaxFactor" -> tauMaxFactor,
            "BinsPerDecade" -> binsPerDecade,
            "WeightPower" -> weightPower,
            "WeightPowerResistance" -> weightPowerResistance,
            "MinFreqRsFit" -> minFreqRsFit,
            "MaxFreqRsFit" -> maxFreqRsFit
          ],
          $Failed
        ],
        {NMinimize::dinfeas}
      ];

      ok = AssociationQ[spec] &&
           KeyExistsQ[spec, "Tau"] && KeyExistsQ[spec, "g"] &&
           VectorQ[spec["Tau"], NumericQ] && VectorQ[spec["g"], NumericQ] &&
           Length[spec["Tau"]] == Length[spec["g"]];

      If[!ok,
        msg = If[spec === $Failed,
          "DCTSpectrum returned $Failed",
          "DCTSpectrum did not return a valid Association"
        ];
      ];

      <|"File" -> file, "Spec" -> spec, "OK" -> ok, "Message" -> msg|>
    ],
    {f, txtFiles}
  ],
  Column[{
    Style["DCT batch progress", 14, Bold],
    Row[{
      ProgressIndicator[i/nFiles, {0, 1}],
      "  ",
      NumberForm[100. i/nFiles, {3, 1}], "%   (", i, "/", nFiles, ")"
    }],
    Row[{"Current file:  ", Style[currentFile, 12]}]
  }]
];

nOK = Count[results[[All, "OK"]], True];
Print["Succeeded on ", nOK, " / ", Length[txtFiles], " files."];

(* ============================================================ *)
(* 3) PLOT RESULTS (rainbow by file order)                       *)
(* ============================================================ *)

goodSpecs = Select[results, #OK === True &];

If[Length[goodSpecs] == 0,
  Print["No successful fits to plot."];
  Abort[];
];

traces = (Transpose[{#Spec["Tau"], #Spec["g"]}] &) /@ goodSpecs;
labels = FileNameTake /@ (goodSpecs[[All, "File"]]);

(* rainbow colors in the same order as traces *)
colors = ColorData["Rainbow"] /@ Rescale[Range[Length[traces]]];

Show[
  ListLinePlot[
    traces,
    PlotStyle -> colors,
    ScalingFunctions -> {"Log10", None},
    Frame -> True,
    FrameLabel -> {"\[Tau] (s)", "g(\[Tau]) (F/decade)"},
    PlotRange -> {Automatic, All},
    PlotLegends -> Placed[labels, Right],
    ImageSize -> 700
  ]
];

(* ============================================================ *)
(* OPTIONAL: PRINT FAILURES                                      *)
(* ============================================================ *)

failures = Select[results, #OK =!= True &];
If[Length[failures] > 0,
  Print["\nFailures:"];
  Do[
    Print["  ", FileNameTake[fail["File"]], " : ", fail["Message"]],
    {fail, failures}
  ];
];


(* ============================================================ *)
(* Plot in electron-transfer-rate space: k = 1/tau (s^-1)        *)
(* ============================================================ *)

tracesK = (Transpose[{1/#Spec["Tau"], #Spec["g"]}] &) /@ goodSpecs;

Show[
  ListLinePlot[
    tracesK,
    PlotStyle -> colors,
    ScalingFunctions -> {"Log10", None},
    Frame -> True,
    FrameLabel -> {"k (s^-1)", "g(k) (F/decade)"},
    PlotRange -> {Automatic, {0, 4 10^-6}},
    PlotLegends -> Placed[labels, Right],
    ImageSize -> 700
  ]
]


(* ::InheritFromParent:: *)
(**)


(#Spec["CPEAlpha"] &) /@ goodSpecs
(#Spec["Y0"] &) /@ goodSpecs
(#Spec["C0"] &) /@ goodSpecs
(#Spec["Rs"] &) /@ goodSpecs
Min[(1/#Spec["Tau"] &) /@ goodSpecs]
Max[(1/#Spec["Tau"] &) /@ goodSpecs]
Min[(1/#Spec["Tau"] &) /@ goodSpecs]/(2 \[Pi])
Max[(1/#Spec["Tau"] &) /@ goodSpecs]/(2 \[Pi])


(* ============================================================ *)
(* 2D heat map in k-space:                                      *)
(*   x = experiment number (1..nExp)                            *)
(*   y = k (s^-1)                                               *)
(*   color = g (F/decade)                                       *)
(* ============================================================ *)

nExp = Length[goodSpecs];
If[nExp == 0, Print["No goodSpecs available."]; Abort[]];

(* ---------- Build a common k-grid across all experiments ---------- *)
kMin = Min[Flatten[(1/#Spec["Tau"]) & /@ goodSpecs]];
kMax = Max[Flatten[(1/#Spec["Tau"]) & /@ goodSpecs]];

nK = 100;  (* resolution in k-direction; increase to 400 if you want smoother *)
kGrid = 10.^Subdivide[Log10[kMin], Log10[kMax], nK - 1];

(* ---------- Interpolate each experiment onto kGrid ---------- *)
gGrid = Table[
   Module[{spec, kData, gData, ord, interp},
     spec = goodSpecs[[i, "Spec"]];
     kData = 1/spec["Tau"];
     gData = spec["g"];

     (* Ensure monotone increasing k for interpolation *)
     ord = Ordering[kData];
     kData = kData[[ord]];
     gData = gData[[ord]];

     (* Piecewise-linear interpolation in log-k space is usually most stable *)
     interp = Interpolation[
       Transpose[{Log10[kData], gData}],
       InterpolationOrder -> 1
     ];

     (* Evaluate on the common grid (with clipping at ends) *)
     interp /@ Log10[kGrid]
   ],
   {i, 1, nExp}
];


expTimesHr = (#[["Spec", "FinishTimeS"]]/3600) & /@ goodSpecs;


(* gGrid is nExp x nK; build points (x,y)->g for ListDensityPlot *)
pts = Flatten[
   Table[
     {expTimesHr[[i]], kGrid[[j]], gGrid[[i, j]]},
     {i, 1, nExp}, {j, 1, nK}
   ],
   1
];

(* ---------- Styling knobs ---------- *)
lblSize = 18;
tickSize = 14;
titleSize = 16;

ListDensityPlot[
  pts,

  (* --- make the whole graphic white --- *)
  Background -> White,
  PlotRangePadding -> Scaled[0.02],

  Frame -> True,
  FrameStyle -> Directive[Black, AbsoluteThickness[1.2]],
  FrameLabel -> {
    Style["Time (hours)", lblSize, Black],
    Style["k (s^-1)",     lblSize, Black]
  },

  (* Log axis on k *)
  ScalingFunctions -> {None, "Log10"},

  (* readable ticks and labels *)
  BaseStyle -> {FontFamily -> "Arial", tickSize, Black},
  LabelStyle -> Directive[lblSize, Black],
  FrameTicksStyle -> Directive[tickSize, Black],

  (* your ticks, but with black tick labels *)
  FrameTicks -> {
    {Automatic, None},  (* left/right ticks for y (k) *)
    {Automatic, None}   (* bottom/top ticks for x (time) *)
  },

  PlotLegends -> Placed[
    BarLegend[Automatic, LabelStyle -> Directive[14, Black],
      LegendMarkerSize -> 220
    ],
    Right
  ],

  InterpolationOrder -> 0,
  ImageSize -> 900,

  PlotLabel -> Style[
    Row[{
      "g(k) heat map   nExp=", nExp,
      "   k-range=[", ScientificForm[kMin, 3], ", ", ScientificForm[kMax, 3], "] s^-1"
    }],
    titleSize, Black
  ]
]



(* ============================================================ *)
(* Rs vs time                                                   *)
(* ============================================================ *)

rsVsTime = Table[
	{
		goodSpecs[[i, "Spec", "FinishTimeS"]]/3600.,
		goodSpecs[[i, "Spec", "Rs"]]
	},
	{i, Length[goodSpecs]}
];

(* ============================================================ *)
(* C0 vs time                                                   *)
(* ============================================================ *)

c0VsTime = Table[
	{
		goodSpecs[[i, "Spec", "FinishTimeS"]]/3600.,
		goodSpecs[[i, "Spec", "C0"]]
	},
	{i, Length[goodSpecs]}
];

GraphicsGrid[
	{
		{
			ListLinePlot[
				rsVsTime,
				Background -> White,
				Frame -> True,
				FrameStyle -> Directive[Black, AbsoluteThickness[1.2]],
				FrameLabel -> {
					Style["Time (hours)", 18, Black],
					Style["Solution resistance Rs (\[CapitalOmega])", 18, Black]
				},
				BaseStyle -> {FontFamily -> "Arial", 14, Black},
				LabelStyle -> Directive[18, Black],
				FrameTicksStyle -> Directive[14, Black],
				PlotStyle -> Directive[Thick],
				PlotMarkers -> {Automatic, 8},
				PlotRange -> {0, Automatic},
				PlotRangePadding -> Scaled[0.02],
				PlotLabel -> Style["Rs vs time", 16, Black]
			],

			ListLinePlot[
				c0VsTime,
				Background -> White,
				Frame -> True,
				FrameStyle -> Directive[Black, AbsoluteThickness[1.2]],
				FrameLabel -> {
					Style["Time (hours)", 18, Black],
					Style["C0 (F)", 18, Black]
				},
				BaseStyle -> {FontFamily -> "Arial", 14, Black},
				LabelStyle -> Directive[18, Black],
				FrameTicksStyle -> Directive[14, Black],
				PlotStyle -> Directive[Thick],
				PlotMarkers -> {Automatic, 8},
				PlotRange -> {0, Automatic},
				PlotRangePadding -> Scaled[0.02],
				PlotLabel -> Style["C0 vs time", 16, Black]
			]
		}
	}, ImageSize -> 900
]


(* ============================================================ *)
(* FAST + MORE ACCURATE: peak picking with sub-grid refinement  *)
(* ============================================================ *)

fitKMin = 10.;
fitKMax = 1000.;

smoothRadius = 1;      (* 0,1,2 ... *)
minSepDec = 0.20;      (* minimum separation in log10(k) *)
minPeakFrac = 0.01;    (* peak must exceed this fraction of max signal *)

(* ---------- helper: simple smoothing ---------- *)
smooth1D[v_List, r_Integer?NonNegative] := Module[{ker},
	If[r == 0, Return[N[v]]];
	ker = ConstantArray[1./(2 r + 1), 2 r + 1];
	ListCorrelate[ker, N[v], {r + 1, -(r + 1)}, 0]
];

(* ---------- helper: local maxima indices ---------- *)
localMaxIndices[y_List] := Select[
	Range[2, Length[y] - 1],
	y[[#]] >= y[[# - 1]] && y[[#]] >= y[[# + 1]] &
];

(* ---------- helper: quadratic peak refinement on 3 points ---------- *)
refinePeakQuadratic[x_List, y_List, idx_Integer] := Module[
	{x1, x2, x3, y1, y2, y3, denom, a, b, xPeak, yPeak},

	If[idx <= 1 || idx >= Length[x],
		Return[<|"xPeak" -> x[[idx]], "yPeak" -> y[[idx]], "Refined" -> False|>]
	];

	x1 = N[x[[idx - 1]]]; x2 = N[x[[idx]]]; x3 = N[x[[idx + 1]]];
	y1 = N[y[[idx - 1]]]; y2 = N[y[[idx]]]; y3 = N[y[[idx + 1]]];

	denom = (x1 - x2) (x1 - x3) (x2 - x3);
	If[denom == 0,
		Return[<|"xPeak" -> x[[idx]], "yPeak" -> y[[idx]], "Refined" -> False|>]
	];

	a = (x3 (y2 - y1) + x2 (y1 - y3) + x1 (y3 - y2))/denom;
	b = (x3^2 (y1 - y2) + x2^2 (y3 - y1) + x1^2 (y2 - y3))/denom;

	If[!NumericQ[a] || !NumericQ[b] || a >= 0,
		Return[<|"xPeak" -> x[[idx]], "yPeak" -> y[[idx]], "Refined" -> False|>]
	];

	xPeak = -b/(2. a);

	If[xPeak < Min[x1, x3] || xPeak > Max[x1, x3],
		Return[<|"xPeak" -> x[[idx]], "yPeak" -> y[[idx]], "Refined" -> False|>]
	];

	yPeak = a xPeak^2 + b xPeak +
		(y1 - a x1^2 - b x1);

	<|"xPeak" -> xPeak, "yPeak" -> yPeak, "Refined" -> True|>
];

(* ---------- helper: interpolate x where y crosses target ---------- *)
interpCrossingX[x1_, y1_, x2_, y2_, yTarget_] := Module[{t},
	If[!And @@ (NumericQ /@ {x1, y1, x2, y2, yTarget}), Return[Missing["NonNumeric"]]];
	If[y2 == y1, Return[Missing["FlatSegment"]]];
	t = (yTarget - y1)/(y2 - y1);
	x1 + t (x2 - x1)
];

(* ---------- helper: FWHM from interpolated half-height crossings ---------- *)
estimateFWHMFromHalfHeightRefined[x_List, y_List, idx_Integer, yPeakOverride_: Automatic] := Module[
	{ypeak, half, iL, iR, xL, xR},

	If[idx < 1 || idx > Length[x], Return[Missing["BadIndex"]]];

	ypeak = If[yPeakOverride === Automatic, y[[idx]], yPeakOverride];
	If[!NumericQ[ypeak] || ypeak <= 0, Return[Missing["BadPeak"]]];

	half = ypeak/2.0;

	iL = idx;
	While[iL > 1 && y[[iL]] > half, iL--];

	iR = idx;
	While[iR < Length[y] && y[[iR]] > half, iR++];

	If[iL == 1 || iR == Length[y], Return[Missing["NoHalfHeightCrossing"]]];

	xL = interpCrossingX[x[[iL]], y[[iL]], x[[iL + 1]], y[[iL + 1]], half];
	xR = interpCrossingX[x[[iR - 1]], y[[iR - 1]], x[[iR]], y[[iR]], half];

	If[!NumericQ[xL] || !NumericQ[xR], Return[Missing["BadCrossing"]]];

	xR - xL
];

(* ---------- helper: choose two separated strongest peaks ---------- *)
pickTwoPeaks[xDec_List, y_List, minSep_?NumericQ, minFrac_?NumericQ] := Module[
	{cand, yMax, goodCand, idx1, remaining, idx2},
	cand = localMaxIndices[y];
	If[cand === {}, Return[{}]];

	yMax = Max[y];
	goodCand = Select[cand, y[[#]] >= minFrac yMax &];
	If[goodCand === {}, Return[{}]];

	idx1 = First @ Ordering[y[[goodCand]], -1];
	idx1 = goodCand[[idx1]];

	remaining = Select[goodCand, Abs[xDec[[#]] - xDec[[idx1]]] >= minSep &];
	If[remaining === {}, Return[{idx1}]];

	idx2 = First @ Ordering[y[[remaining]], -1];
	idx2 = remaining[[idx2]];

	SortBy[{idx1, idx2}, xDec[[#]] &]
];

(* ---------- analyze one trace ---------- *)
analyzeOneTraceFastRefined[selector_Integer] := Module[
	{
		kUse0, gUse0, xDec, yRaw, ySm, peakIdx, idx1, idx2,
		ref1, ref2, xPeakDec1, xPeakDec2, yPeak1, yPeak2,
		kPeak1, kPeak2, fwhmDec1, fwhmDec2, fwhmK1, fwhmK2
	},

	kUse0 = Pick[kGrid, Map[fitKMin <= # <= fitKMax &, kGrid]];
	gUse0 = Pick[gGrid[[selector, All]], Map[fitKMin <= # <= fitKMax &, kGrid]];

	If[Length[kUse0] < 8,
		Return[<|"OK" -> False, "Selector" -> selector, "Message" -> "Too few points in fit range."|>]
	];

	xDec = Log10[kUse0];
	yRaw = N[gUse0];
	ySm = smooth1D[yRaw, smoothRadius];
	peakIdx = pickTwoPeaks[xDec, ySm, minSepDec, minPeakFrac];

	If[Length[peakIdx] < 2,
		Return[<|"OK" -> False, "Selector" -> selector, "Message" -> "Could not find two separated peaks."|>]
	];

	idx1 = peakIdx[[1]];
	idx2 = peakIdx[[2]];

	ref1 = refinePeakQuadratic[xDec, ySm, idx1];
	ref2 = refinePeakQuadratic[xDec, ySm, idx2];

	xPeakDec1 = ref1["xPeak"];
	xPeakDec2 = ref2["xPeak"];
	yPeak1 = ref1["yPeak"];
	yPeak2 = ref2["yPeak"];

	kPeak1 = 10.^xPeakDec1;
	kPeak2 = 10.^xPeakDec2;

	fwhmDec1 = estimateFWHMFromHalfHeightRefined[xDec, ySm, idx1, yPeak1];
	fwhmDec2 = estimateFWHMFromHalfHeightRefined[xDec, ySm, idx2, yPeak2];

	fwhmK1 = If[NumericQ[fwhmDec1],
		10.^(xPeakDec1 + fwhmDec1/2) - 10.^(xPeakDec1 - fwhmDec1/2),
		Missing["NoFWHM"]
	];

	fwhmK2 = If[NumericQ[fwhmDec2],
		10.^(xPeakDec2 + fwhmDec2/2) - 10.^(xPeakDec2 - fwhmDec2/2),
		Missing["NoFWHM"]
	];

	<|
		"OK" -> True,
		"Selector" -> selector,
		"TimeHr" -> expTimesHr[[selector]],

		"kGridUsed" -> kUse0,
		"xDecUsed" -> xDec,
		"yRaw" -> yRaw,
		"ySm" -> ySm,

		"idx1" -> idx1,
		"idx2" -> idx2,

		"xPeakDec1" -> xPeakDec1,
		"xPeakDec2" -> xPeakDec2,
		"yPeak1" -> yPeak1,
		"yPeak2" -> yPeak2,

		"kPeak1" -> kPeak1,
		"kPeak2" -> kPeak2,

		"FWHMDec1" -> fwhmDec1,
		"FWHMDec2" -> fwhmDec2,
		"FWHMK1" -> fwhmK1,
		"FWHMK2" -> fwhmK2
	|>
];

(* ---------- run across all traces with progress ---------- *)
nTraces = Length[gGrid];
fitCounter = 0;
currentSelector = 0;

allPeakPicks = ConstantArray[Missing["NotComputed"], nTraces];

Monitor[
	Do[
		currentSelector = i;
		fitCounter = i;
		allPeakPicks[[i]] = analyzeOneTraceFastRefined[i];
	,
		{i, 1, nTraces}
	],
	Column[{
		Style["Fast refined peak-picking progress", 14, Bold],
		Row[{
			ProgressIndicator[N[fitCounter/nTraces], {0, 1}],
			"  ",
			NumberForm[100.0 fitCounter/nTraces, {4, 1}],
			"%   (", fitCounter, "/", nTraces, ")"
		}],
		Row[{"Current selector: ", currentSelector}]
	}]
];

goodPeakPicks = Select[allPeakPicks, AssociationQ[#] && TrueQ[#["OK"]] &];
badPeakPicks  = Select[allPeakPicks, AssociationQ[#] && !TrueQ[#["OK"]] &];

Print["Succeeded on ", Length[goodPeakPicks], " / ", nTraces, " traces."];

If[Length[badPeakPicks] > 0,
	Print["Failures:"];
	Scan[
		(Print["  selector ", #["Selector"], ": ", #["Message"]]) &,
		badPeakPicks
	];
];

If[Length[goodPeakPicks] == 0,
	Print["No successful peak picks."];
	Abort[];
];

(* ---------- build time series ---------- *)
peak1AmpVsTime = Table[
	{goodPeakPicks[[i, "TimeHr"]], goodPeakPicks[[i, "yPeak1"]]},
	{i, Length[goodPeakPicks]}
];

peak2AmpVsTime = Table[
	{goodPeakPicks[[i, "TimeHr"]], goodPeakPicks[[i, "yPeak2"]]},
	{i, Length[goodPeakPicks]}
];


peak1VsTime = Table[
	{goodPeakPicks[[i, "TimeHr"]], goodPeakPicks[[i, "kPeak1"]]},
	{i, Length[goodPeakPicks]}
];

peak2VsTime = Table[
	{goodPeakPicks[[i, "TimeHr"]], goodPeakPicks[[i, "kPeak2"]]},
	{i, Length[goodPeakPicks]}
];

fwhm1VsTime = Select[
	Table[
		{goodPeakPicks[[i, "TimeHr"]], goodPeakPicks[[i, "FWHMK1"]]},
		{i, Length[goodPeakPicks]}
	],
	NumericQ[#[[2]]] && #[[2]] > 0 &
];

fwhm2VsTime = Select[
	Table[
		{goodPeakPicks[[i, "TimeHr"]], goodPeakPicks[[i, "FWHMK2"]]},
		{i, Length[goodPeakPicks]}
	],
	NumericQ[#[[2]]] && #[[2]] > 0 &
];

(* ---------- plots: time evolution ---------- *)
peakAmplitudesPlot =
	ListPlot[
		{peak1AmpVsTime, peak2AmpVsTime},
		Joined -> True,
		PlotStyle -> {Directive[Blue, Thick], Directive[Red, Thick]},
		PlotMarkers -> {{Automatic, 7}, {Automatic, 7}},
		Frame -> True,
		Axes -> False,
		Background -> White,
		ImageSize -> 500,
		FrameLabel -> {"Time (hours)", "Peak Amplitude (F)"},
		PlotRange -> All,
		PlotLegends -> Placed[{"Peak 1", "Peak 2"}, Right],
		PlotLabel -> "Peak amplitude vs time"
	];


peakCentersPlot =
	ListPlot[
		{peak1VsTime, peak2VsTime},
		Joined -> True,
		PlotStyle -> {Directive[Blue, Thick], Directive[Red, Thick]},
		PlotMarkers -> {{Automatic, 7}, {Automatic, 7}},
		Frame -> True,
		Axes -> False,
		Background -> White,
		ImageSize -> 500,
		FrameLabel -> {"Time (hours)", "Peak center k (s^-1)"},
		PlotRange -> {Automatic,{0,300}},
		PlotLegends -> Placed[{"Peak 1", "Peak 2"}, Right],
		PlotLabel -> "Peak k vs time"
	];

fwhmPlot =
	ListPlot[
		{fwhm1VsTime, fwhm2VsTime},
		Joined -> True,
		PlotStyle -> {Directive[Blue, Thick], Directive[Red, Thick]},
		PlotMarkers -> {{Automatic, 7}, {Automatic, 7}},
		Frame -> True,
		Axes -> False,
		Background -> White,
		ImageSize -> 500,
		FrameLabel -> {"Time (hours)", "Estimated FWHM in k (s^-1)"},
		PlotRange -> All,
		PlotLegends -> Placed[{"Peak 1", "Peak 2"}, Right],
		PlotLabel -> "Peak FWHM vs time"
	];
peakAmplitudesPlot
peakCentersPlot
fwhmPlot



exampleSelector = 1;   (* choose one trace to visualize *)

exampleResult = analyzeOneTraceFastRefined[exampleSelector];

If[!TrueQ@exampleResult["OK"],
	Print["Example selector failed: ", exampleResult["Message"]],
	
	Module[
		{
			kUse, xDec, yRaw,
			xpk1, xpk2, ypk1, ypk2, fwhm1, fwhm2,
			sigma1, sigma2, fitFun, fitVals,
			rawPts, fitPts, valid
		},
		
		kUse = exampleResult["kGridUsed"];
		xDec = exampleResult["xDecUsed"];
		yRaw = exampleResult["yRaw"];
		
		xpk1 = exampleResult["xPeakDec1"];
		xpk2 = exampleResult["xPeakDec2"];
		ypk1 = exampleResult["yPeak1"];
		ypk2 = exampleResult["yPeak2"];
		
		fwhm1 = exampleResult["FWHMDec1"];
		fwhm2 = exampleResult["FWHMDec2"];
		
		sigma1 = If[NumericQ[fwhm1] && fwhm1 > 0,
			fwhm1/(2 Sqrt[2 Log[2]]),
			Missing["BadSigma1"]
		];
		
		sigma2 = If[NumericQ[fwhm2] && fwhm2 > 0,
			fwhm2/(2 Sqrt[2 Log[2]]),
			Missing["BadSigma2"]
		];
		
		valid = And[
			VectorQ[kUse, NumericQ],
			VectorQ[yRaw, NumericQ],
			NumericQ[xpk1], NumericQ[xpk2],
			NumericQ[ypk1], NumericQ[ypk2],
			NumericQ[sigma1], NumericQ[sigma2]
		];
		
		If[!valid,
			Print["Could not build example plot because some fit quantities are nonnumeric."],
			
			fitFun[x_] := 
				ypk1 Exp[-(x - xpk1)^2/(2 sigma1^2)] +
				ypk2 Exp[-(x - xpk2)^2/(2 sigma2^2)];
			
			fitVals = fitFun /@ xDec;
			
			rawPts = Transpose[{kUse, yRaw}];
			fitPts = Transpose[{kUse, fitVals}];
			
			Print[
				ListPlot[
					{rawPts, fitPts},
					Joined -> {False, True},
					PlotStyle -> {
						Directive[Black, PointSize[0.015]],
						Directive[Red, Thick]
					},
					Frame -> True,
					Axes -> False,
					Background -> White,
					ImageSize -> 900,
					ScalingFunctions -> {"Log10", None},
					FrameLabel -> {"k (s^-1)", "g(k)"},
					PlotRange -> All,
					PlotLegends -> Placed[{"Raw data", "Reconstructed fit"}, Right],
					PlotLabel -> Row[{
						"Selected trace and reconstructed fit, selector = ",
						exampleSelector,
						", time = ",
						NumberForm[exampleResult["TimeHr"], {6, 2}],
						" hr"
					}]
				]
			];
		];
	]
];


(* ============================================================ *)
(* Rainbow Bode overlays: RAW data for ALL files                *)
(*   Col 1: RAW                                                 *)
(*   Col 2: FIT                                                 *)
(* ============================================================ *)

ImportEISRaw[file_String] := Module[
  {lines, rows, nums, good, f, zre, zimNeg, z, colF, colZre, colZim},

  lines = Import[file, "Lines"];
  lines = Select[lines, StringTrim[#] =!= "" &];

  rows = StringSplit[#, ";"] & /@ lines;
  nums = (Quiet@Check[ToExpression[StringTrim[#]], Missing["Bad"]] &) /@ rows;

  (* keep rows whose first 4 entries are numeric *)
  good = Select[nums, Length[#] >= 4 && And @@ (NumericQ /@ #[[1 ;; 4]]) &];
  If[good === {}, Return[$Failed]];

  (* Detect whether col1 is Index (integer-ish) by checking monotone small integers *)
  (* If col1 looks like 1,2,3,... then Frequency is column 2; else Frequency is column 1 *)
  If[
    VectorQ[good[[All, 1]], (NumericQ[#] && # >= 0 && # <= 10^6) &] &&
      (* heuristic: many index columns start at 1 and increment by 1 *)
      (Max[Abs[Differences[good[[All, 1]]]]] <= 2),
    colF = 2; colZre = 3; colZim = 4,
    colF = 1; colZre = 2; colZim = 3
  ];

  f      = good[[All, colF]];
  zre    = good[[All, colZre]];
  zimNeg = good[[All, colZim]];
  z      = zre + I*(-zimNeg);  (* file stores -Z'' *)

  (* Sort by frequency ascending for consistent plotting *)
  With[{ord = Ordering[f]},
    <|"FreqHz" -> f[[ord]], "Z" -> z[[ord]]|>
  ]
];



(* ============================================================ *)
(* Impedance Bode overlays                                      *)
(*   Col 1: RAW impedance                                       *)
(*   Col 2: FIT impedance                                       *)
(*   Col 3: ABS residual error                                  *)
(*   Row 1: |Z|                                                 *)
(*   Row 2: Phase(Z)                                            *)
(*   Row 3: Re(Z)                                               *)
(*   Row 4: -Im(Z)                                              *)
(* accommodates useConstantPhaseElementFlag                     *)
(* ============================================================ *)

If[Length[goodSpecs] == 0,
	Print["No goodSpecs to plot."];
	Abort[];
];

labels = FileNameTake /@ (goodSpecs[[All, "File"]]);
files  = goodSpecs[[All, "File"]];

(* ---------- RAW impedance ---------- *)
rawZResults =
	Table[
		Module[
			{
				file, raw, fRaw, zRaw,
				magTrace, phaseTrace, reTrace, imTrace, ok, msg
			},

			file = files[[k]];
			raw  = Quiet @ Check[ImportEISRaw[file], $Failed];

			If[
				raw === $Failed || !AssociationQ[raw] || !KeyExistsQ[raw, "FreqHz"] || !KeyExistsQ[raw, "Z"],
				Return[
					<|
						"File" -> file,
						"Label" -> labels[[k]],
						"OK" -> False,
						"Message" -> "raw parse failed"
					|>
				]
			];

			fRaw = raw["FreqHz"];
			zRaw = raw["Z"];

			magTrace   = Transpose[{fRaw, Abs[zRaw]}];
			phaseTrace = Transpose[{fRaw, (180./Pi) * Arg[zRaw]}];
			reTrace    = Transpose[{fRaw, Re[zRaw]}];
			imTrace    = Transpose[{fRaw, -Im[zRaw]}];

			ok = VectorQ[magTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &] &&
				 VectorQ[phaseTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &] &&
				 VectorQ[reTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &] &&
				 VectorQ[imTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &];

			msg = If[ok, "", "raw impedance trace non-numeric"];

			<|
				"File" -> file,
				"Label" -> labels[[k]],
				"FreqHz" -> fRaw,
				"ZRaw" -> zRaw,
				"OK" -> ok,
				"MagTrace" -> magTrace,
				"PhaseTrace" -> phaseTrace,
				"ReTrace" -> reTrace,
				"ImTrace" -> imTrace,
				"Message" -> msg
			|>
		],
		{k, Length[files]}
	];

goodRawZ = Select[rawZResults, TrueQ[#["OK"]] &];
badRawZ  = Select[rawZResults, !TrueQ[#["OK"]] &];

Print["Raw impedance overlay: parsed OK = ", Length[goodRawZ], " / ", Length[rawZResults]];

If[Length[badRawZ] > 0,
	Print["Failures:"];
	Scan[
		(Print["\t", FileNameTake[#["File"]], " : ", #["Message"]]) &,
		badRawZ
	];
];

If[Length[goodRawZ] == 0,
	Print["No raw impedance spectra could be parsed."];
	Abort[];
];

rawZMagTraces   = goodRawZ[[All, "MagTrace"]];
rawZPhaseTraces = goodRawZ[[All, "PhaseTrace"]];
rawZReTraces    = goodRawZ[[All, "ReTrace"]];
rawZImTraces    = goodRawZ[[All, "ImTrace"]];
rawZLabels      = goodRawZ[[All, "Label"]];
rawZColors      = ColorData["Rainbow"] /@ Rescale[Range[Length[goodRawZ]]];

(* ---------- FIT impedance ---------- *)
fullFitZResults =
	Table[
		Module[
			{
				spec, raw, tau, g, c0, y0, alphaCPE, useCPE, rs,
				fRaw, \[Omega], \[CapitalDelta]log, Kmat, YintFit, Zfit,
				magTrace, phaseTrace, reTrace, imTrace, file, ok, msg
			},

			spec = goodSpecs[[k, "Spec"]];
			file = goodSpecs[[k, "File"]];
			raw  = Quiet @ Check[ImportEISRaw[file], $Failed];

			If[
				raw === $Failed || !AssociationQ[raw] || !KeyExistsQ[raw, "FreqHz"],
				Return[
					<|
						"File" -> file,
						"Label" -> labels[[k]],
						"OK" -> False,
						"Message" -> "raw parse failed or missing FreqHz"
					|>
				]
			];

			tau = spec["Tau"];
			g   = spec["g"];
			c0  = Lookup[spec, "C0", Missing["NotFound"]];
			y0  = Lookup[spec, "Y0", Missing["NotFound"]];
			alphaCPE = Lookup[spec, "CPEAlpha", Missing["NotFound"]];
			useCPE = TrueQ[Lookup[spec, "useConstantPhaseElementFlag", False]];
			rs = spec["Rs"];

			fRaw = raw["FreqHz"];
			\[Omega] = 2 Pi fRaw;

			\[CapitalDelta]log = Mean[Differences[Log10[tau]]];

			Kmat = Table[
				(I*\[Omega][[j]])/(1 + I*\[Omega][[j]]*tau[[m]]),
				{j, Length[\[Omega]]}, {m, Length[tau]}
			];

			YintFit =
				If[useCPE,
					If[!NumericQ[y0] || !NumericQ[alphaCPE],
						$Failed,
						y0*(I*\[Omega])^alphaCPE + Kmat . (g * \[CapitalDelta]log)
					],
					If[!NumericQ[c0],
						$Failed,
						(I*\[Omega])*c0 + Kmat . (g * \[CapitalDelta]log)
					]
				];

			If[YintFit === $Failed,
				Return[
					<|
						"File" -> file,
						"Label" -> labels[[k]],
						"OK" -> False,
						"Message" -> "missing fit parameters for impedance reconstruction"
					|>
				]
			];

			If[Min[Abs[YintFit]] < 10^-30,
				Return[
					<|
						"File" -> file,
						"Label" -> labels[[k]],
						"OK" -> False,
						"Message" -> "YintFit too small"
					|>
				]
			];

			Zfit = rs + 1/YintFit;

			magTrace   = Transpose[{fRaw, Abs[Zfit]}];
			phaseTrace = Transpose[{fRaw, (180./Pi) * Arg[Zfit]}];
			reTrace    = Transpose[{fRaw, Re[Zfit]}];
			imTrace    = Transpose[{fRaw, -Im[Zfit]}];

			ok = VectorQ[magTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &] &&
				 VectorQ[phaseTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &] &&
				 VectorQ[reTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &] &&
				 VectorQ[imTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &];

			msg = If[ok, "", "fit impedance trace non-numeric"];

			<|
				"File" -> file,
				"Label" -> labels[[k]],
				"FreqHz" -> fRaw,
				"ZFit" -> Zfit,
				"MagTrace" -> magTrace,
				"PhaseTrace" -> phaseTrace,
				"ReTrace" -> reTrace,
				"ImTrace" -> imTrace,
				"OK" -> ok,
				"Message" -> msg
			|>
		],
		{k, Length[goodSpecs]}
	];

goodFullFitZ = Select[fullFitZResults, TrueQ[#["OK"]] &];
badFullFitZ  = Select[fullFitZResults, !TrueQ[#["OK"]] &];

If[Length[badFullFitZ] > 0,
	Print["Failures:"];
	Scan[
		(Print["\t", #["File"], " : ", #["Message"]]) &,
		badFullFitZ
	];
];

If[Length[goodFullFitZ] == 0,
	Print["No fit impedance spectra could be computed."];
	Abort[];
];

fullFitZMagTraces   = goodFullFitZ[[All, "MagTrace"]];
fullFitZPhaseTraces = goodFullFitZ[[All, "PhaseTrace"]];
fullFitZReTraces    = goodFullFitZ[[All, "ReTrace"]];
fullFitZImTraces    = goodFullFitZ[[All, "ImTrace"]];
fitZColors          = ColorData["Rainbow"] /@ Rescale[Range[Length[goodFullFitZ]]];

(* ---------- Residual traces ---------- *)
residualZResults =
	Table[
		Module[
			{
				rawAssoc, fitAssoc, fRaw, zRaw, zFit,
				magResTrace, phaseResTrace, reResTrace, imResTrace
			},

			rawAssoc = rawZResults[[k]];
			fitAssoc = fullFitZResults[[k]];

			If[!TrueQ[rawAssoc["OK"]] || !TrueQ[fitAssoc["OK"]],
				Return[
					<|
						"OK" -> False,
						"Message" -> "raw/fit pair unavailable"
					|>
				]
			];

			fRaw = rawAssoc["FreqHz"];
			zRaw = rawAssoc["ZRaw"];
			zFit = fitAssoc["ZFit"];

			If[Length[fRaw] =!= Length[zRaw] || Length[zRaw] =!= Length[zFit],
				Return[
					<|
						"OK" -> False,
						"Message" -> "raw and fit lengths mismatch"
					|>
				]
			];

			magResTrace   = Transpose[{fRaw, Abs[(zRaw - zFit)/zRaw]}];
			phaseResTrace = Transpose[{fRaw, Abs[(180./Pi) Arg[zRaw] - (180./Pi) Arg[zFit]]}];
			reResTrace    = Transpose[{fRaw, Abs[Re[zRaw] - Re[zFit]]}];
			imResTrace    = Transpose[{fRaw, Abs[(-Im[zRaw]) - (-Im[zFit])]}];

			<|
				"OK" -> True,
				"MagResTrace" -> magResTrace,
				"PhaseResTrace" -> phaseResTrace,
				"ReResTrace" -> reResTrace,
				"ImResTrace" -> imResTrace
			|>
		],
		{k, Length[files]}
	];

goodResidualZ = Select[residualZResults, TrueQ[#["OK"]] &];

If[Length[goodResidualZ] == 0,
	Print["No residual impedance traces could be computed."];
	Abort[];
];

residualZMagTraces   = goodResidualZ[[All, "MagResTrace"]];
residualZPhaseTraces = goodResidualZ[[All, "PhaseResTrace"]];
residualZReTraces    = goodResidualZ[[All, "ReResTrace"]];
residualZImTraces    = goodResidualZ[[All, "ImResTrace"]];
resZColors           = ColorData["Rainbow"] /@ Rescale[Range[Length[goodResidualZ]]];

(* ---------- Shared y-axis ranges ---------- *)

allMagY = Select[
	Join[
		Flatten[rawZMagTraces[[All, All, 2]]],
		Flatten[fullFitZMagTraces[[All, All, 2]]]
	],
	NumericQ[#] && # > 0 &
];

allPhaseY = Select[
	Join[
		Flatten[rawZPhaseTraces[[All, All, 2]]],
		Flatten[fullFitZPhaseTraces[[All, All, 2]]]
	],
	NumericQ
];

allReY = Select[
	Join[
		Flatten[rawZReTraces[[All, All, 2]]],
		Flatten[fullFitZReTraces[[All, All, 2]]]
	],
	NumericQ
];

allImY = Select[
	Join[
		Flatten[rawZImTraces[[All, All, 2]]],
		Flatten[fullFitZImTraces[[All, All, 2]]]
	],
	NumericQ
];

allMagResY = Select[Flatten[residualZMagTraces[[All, All, 2]]], NumericQ[#] && # > 0 &];
allPhaseResY = Select[Flatten[residualZPhaseTraces[[All, All, 2]]], NumericQ];
allReResY = Select[Flatten[residualZReTraces[[All, All, 2]]], NumericQ];
allImResY = Select[Flatten[residualZImTraces[[All, All, 2]]], NumericQ];

If[
	allMagY === {} || allPhaseY === {} || allReY === {} || allImY === {} ||
	allMagResY === {} || allPhaseResY === {} || allReResY === {} || allImResY === {},
	Print["Could not determine shared axis ranges."];
	Abort[];
];

magYRange      = {Min[allMagY], Max[allMagY]};
phaseYRange    = {Min[allPhaseY], Max[allPhaseY]};
reYRange       = {Min[allReY], Max[allReY]};
imYRange       = {Min[allImY], Max[allImY]};
magResYRange   = {Min[allMagResY], Max[allMagResY]};
phaseResYRange = {Min[allPhaseResY], Max[allPhaseResY]};
reResYRange    = {Min[allReResY], Max[allReResY]};
imResYRange    = {Min[allImResY], Max[allImResY]};

Print["Shared |Z| y-range = ", magYRange];
Print["Shared phase y-range = ", phaseYRange];
Print["Shared Re(Z) y-range = ", reYRange];
Print["Shared -Im(Z) y-range = ", imYRange];

(* ---------- Plots ---------- *)

rawZMagPlot =
	ListLogLogPlot[
		rawZMagTraces,
		PlotStyle -> rawZColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "|Zraw| (\[CapitalOmega])"},
		PlotRange -> {All, magYRange},
		ImageSize -> 500,
		Background -> White,
		PlotLabel -> Row[{"RAW |Z|   n=", Length[goodRawZ]}]
	];

fitZMagPlot =
	ListLogLogPlot[
		fullFitZMagTraces,
		PlotStyle -> fitZColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "|Zfit| (\[CapitalOmega])"},
		PlotRange -> {All, magYRange},
		ImageSize -> 500,
		Background -> White,
		PlotLabel -> Row[{"FIT |Z|   n=", Length[goodFullFitZ]}]
	];

resZMagPlot =
	ListLogLogPlot[
		residualZMagTraces,
		PlotStyle -> resZColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "|(Zraw - Zfit)/Zraw|"},
		PlotRange -> {All, magResYRange},
		ImageSize -> 500,
		Background -> White,
		PlotLabel -> "ABS residual error |(Zraw - Zfit)/Zraw|"
	];

rawZPhasePlot =
	ListLogLinearPlot[
		rawZPhaseTraces,
		PlotStyle -> rawZColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "Phase(Zraw) (deg)"},
		PlotRange -> {All, phaseYRange},
		ImageSize -> 500,
		Background -> White,
		PlotLabel -> Row[{"RAW phase(Z)   n=", Length[goodRawZ]}]
	];

fitZPhasePlot =
	ListLogLinearPlot[
		fullFitZPhaseTraces,
		PlotStyle -> fitZColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "Phase(Zfit) (deg)"},
		PlotRange -> {All, phaseYRange},
		ImageSize -> 500,
		Background -> White,
		PlotLabel -> Row[{"FIT phase(Z)   n=", Length[goodFullFitZ]}]
	];

resZPhasePlot =
	ListLogLinearPlot[
		residualZPhaseTraces,
		PlotStyle -> resZColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "|Phase raw - fit| (deg)"},
		PlotRange -> {All, phaseResYRange},
		ImageSize -> 500,
		Background -> White,
		PlotLabel -> "ABS phase residual"
	];

rawZRePlot =
	ListLogLogPlot[
		rawZReTraces,
		PlotStyle -> rawZColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "Re(Zraw) (\[CapitalOmega])"},
		PlotRange -> {All, reYRange},
		ImageSize -> 500,
		Background -> White,
		PlotLabel -> "RAW Re(Z)"
	];

fitZRePlot =
	ListLogLogPlot[
		fullFitZReTraces,
		PlotStyle -> fitZColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "Re(Zfit) (\[CapitalOmega])"},
		PlotRange -> {All, reYRange},
		ImageSize -> 500,
		Background -> White,
		PlotLabel -> "FIT Re(Z)"
	];

resZRePlot =
	ListLogLogPlot[
		residualZReTraces,
		PlotStyle -> resZColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "|Re raw - fit| (\[CapitalOmega])"},
		PlotRange -> {All, reResYRange},
		ImageSize -> 500,
		Background -> White,
		PlotLabel -> "ABS Re residual"
	];

rawZImPlot =
	ListLogLogPlot[
		rawZImTraces,
		PlotStyle -> rawZColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "-Im(Zraw) (\[CapitalOmega])"},
		PlotRange -> {All, imYRange},
		ImageSize -> 500,
		Background -> White,
		PlotLabel -> "RAW -Im(Z)"
	];

fitZImPlot =
	ListLogLogPlot[
		fullFitZImTraces,
		PlotStyle -> fitZColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "-Im(Zfit) (\[CapitalOmega])"},
		PlotRange -> {All, imYRange},
		ImageSize -> 500,
		Background -> White,
		PlotLabel -> "FIT -Im(Z)"
	];

resZImPlot =
	ListLogLogPlot[
		residualZImTraces,
		PlotStyle -> resZColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "|(-Im raw) - (-Im fit)| (\[CapitalOmega])"},
		PlotRange -> {All, imResYRange},
		ImageSize -> 500,
		Background -> White,
		PlotLabel -> "ABS -Im residual"
	];

GraphicsGrid[
	{
		{rawZMagPlot, fitZMagPlot, resZMagPlot},
		{rawZPhasePlot, fitZPhasePlot, resZPhasePlot},
		{rawZRePlot, fitZRePlot, resZRePlot},
		{rawZImPlot, fitZImPlot, resZImPlot}
	}, ImageSize->1200, Spacings->{0,0}
]


(* ============================================================ *)
(* Admittance Bode overlays after removing Rs                   *)
(*   Col 1: RAW interface admittance                            *)
(*   Col 2: FIT interface admittance                            *)
(*   Col 3: ABS residual error                                  *)
(*   Row 1: |Yint|                                              *)
(*   Row 2: Phase(Yint)                                         *)
(*   Row 3: Re(Yint)                                            *)
(*   Row 4: -Im(Yint)                                           *)
(* accommodates useConstantPhaseElementFlag                     *)
(* ============================================================ *)

If[Length[goodSpecs] == 0,
	Print["No goodSpecs to plot."];
	Abort[];
];

labels = FileNameTake /@ (goodSpecs[[All, "File"]]);
files  = goodSpecs[[All, "File"]];

(* ---------- RAW admittance after removing Rs ---------- *)
rawYintResults =
	Table[
		Module[
			{
				file, raw, spec, rs, fRaw, zRaw, zIntRaw, yIntRaw,
				magTrace, phaseTrace, reTrace, imTrace, ok, msg
			},

			file = files[[k]];
			raw  = Quiet @ Check[ImportEISRaw[file], $Failed];
			spec = goodSpecs[[k, "Spec"]];

			If[
				raw === $Failed || !AssociationQ[raw] || !KeyExistsQ[raw, "FreqHz"] || !KeyExistsQ[raw, "Z"],
				Return[
					<|
						"File" -> file,
						"Label" -> labels[[k]],
						"OK" -> False,
						"Message" -> "raw parse failed"
					|>
				]
			];

			rs   = spec["Rs"];
			fRaw = raw["FreqHz"];
			zRaw = raw["Z"];

			zIntRaw = zRaw - rs;

			If[Min[Abs[zIntRaw]] < 10^-15,
				Return[
					<|
						"File" -> file,
						"Label" -> labels[[k]],
						"OK" -> False,
						"Message" -> "Zraw - Rs too small"
					|>
				]
			];

			yIntRaw = 1/zIntRaw;

			magTrace   = Transpose[{fRaw, Abs[yIntRaw]}];
			phaseTrace = Transpose[{fRaw, (180./Pi) * Arg[yIntRaw]}];
			reTrace    = Transpose[{fRaw, Re[yIntRaw]}];
			imTrace    = Transpose[{fRaw, -Im[yIntRaw]}];

			ok = VectorQ[magTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &] &&
				 VectorQ[phaseTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &] &&
				 VectorQ[reTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &] &&
				 VectorQ[imTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &];

			msg = If[ok, "", "raw admittance trace non-numeric"];

			<|
				"File" -> file,
				"Label" -> labels[[k]],
				"FreqHz" -> fRaw,
				"YIntRaw" -> yIntRaw,
				"OK" -> ok,
				"MagTrace" -> magTrace,
				"PhaseTrace" -> phaseTrace,
				"ReTrace" -> reTrace,
				"ImTrace" -> imTrace,
				"Message" -> msg
			|>
		],
		{k, Length[files]}
	];

goodRawYint = Select[rawYintResults, TrueQ[#["OK"]] &];
badRawYint  = Select[rawYintResults, !TrueQ[#["OK"]] &];

Print["Raw admittance overlay: parsed OK = ", Length[goodRawYint], " / ", Length[rawYintResults]];

If[Length[badRawYint] > 0,
	Print["Failures:"];
	Scan[
		(Print["\t", FileNameTake[#["File"]], " : ", #["Message"]]) &,
		badRawYint
	];
];

If[Length[goodRawYint] == 0,
	Print["No raw admittance spectra could be parsed."];
	Abort[];
];

rawYintMagTraces   = goodRawYint[[All, "MagTrace"]];
rawYintPhaseTraces = goodRawYint[[All, "PhaseTrace"]];
rawYintReTraces    = goodRawYint[[All, "ReTrace"]];
rawYintImTraces    = goodRawYint[[All, "ImTrace"]];
rawYintLabels      = goodRawYint[[All, "Label"]];
rawYintColors      = ColorData["Rainbow"] /@ Rescale[Range[Length[goodRawYint]]];

(* ---------- FIT admittance after removing Rs ---------- *)
fullFitYintResults =
	Table[
		Module[
			{
				spec, raw, tau, g, c0, y0, alphaCPE, useCPE,
				fRaw, \[Omega], \[CapitalDelta]log, Kmat, YintFit,
				magTrace, phaseTrace, reTrace, imTrace, file, ok, msg
			},

			spec = goodSpecs[[k, "Spec"]];
			file = goodSpecs[[k, "File"]];
			raw  = Quiet @ Check[ImportEISRaw[file], $Failed];

			If[
				raw === $Failed || !AssociationQ[raw] || !KeyExistsQ[raw, "FreqHz"],
				Return[
					<|
						"File" -> file,
						"Label" -> labels[[k]],
						"OK" -> False,
						"Message" -> "raw parse failed or missing FreqHz"
					|>
				]
			];

			tau = spec["Tau"];
			g   = spec["g"];
			c0  = Lookup[spec, "C0", Missing["NotFound"]];
			y0  = Lookup[spec, "Y0", Missing["NotFound"]];
			alphaCPE = Lookup[spec, "CPEAlpha", Missing["NotFound"]];
			useCPE = TrueQ[Lookup[spec, "useConstantPhaseElementFlag", False]];

			fRaw = raw["FreqHz"];
			\[Omega] = 2 Pi fRaw;

			\[CapitalDelta]log = Mean[Differences[Log10[tau]]];

			Kmat = Table[
				(I*\[Omega][[j]])/(1 + I*\[Omega][[j]]*tau[[m]]),
				{j, Length[\[Omega]]}, {m, Length[tau]}
			];

			YintFit =
				If[useCPE,
					If[!NumericQ[y0] || !NumericQ[alphaCPE],
						$Failed,
						y0*(I*\[Omega])^alphaCPE + Kmat . (g * \[CapitalDelta]log)
					],
					If[!NumericQ[c0],
						$Failed,
						(I*\[Omega])*c0 + Kmat . (g * \[CapitalDelta]log)
					]
				];

			If[YintFit === $Failed,
				Return[
					<|
						"File" -> file,
						"Label" -> labels[[k]],
						"OK" -> False,
						"Message" -> "missing fit parameters for admittance reconstruction"
					|>
				]
			];

			magTrace   = Transpose[{fRaw, Abs[YintFit]}];
			phaseTrace = Transpose[{fRaw, (180./Pi) * Arg[YintFit]}];
			reTrace    = Transpose[{fRaw, Re[YintFit]}];
			imTrace    = Transpose[{fRaw, -Im[YintFit]}];

			ok = VectorQ[magTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &] &&
				 VectorQ[phaseTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &] &&
				 VectorQ[reTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &] &&
				 VectorQ[imTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &];

			msg = If[ok, "", "fit admittance trace non-numeric"];

			<|
				"File" -> file,
				"Label" -> labels[[k]],
				"FreqHz" -> fRaw,
				"YIntFit" -> YintFit,
				"MagTrace" -> magTrace,
				"PhaseTrace" -> phaseTrace,
				"ReTrace" -> reTrace,
				"ImTrace" -> imTrace,
				"OK" -> ok,
				"Message" -> msg
			|>
		],
		{k, Length[goodSpecs]}
	];

goodFullFitYint = Select[fullFitYintResults, TrueQ[#["OK"]] &];
badFullFitYint  = Select[fullFitYintResults, !TrueQ[#["OK"]] &];

If[Length[badFullFitYint] > 0,
	Print["Failures:"];
	Scan[
		(Print["\t", #["File"], " : ", #["Message"]]) &,
		badFullFitYint
	];
];

If[Length[goodFullFitYint] == 0,
	Print["No fit admittance spectra could be computed."];
	Abort[];
];

fullFitYintMagTraces   = goodFullFitYint[[All, "MagTrace"]];
fullFitYintPhaseTraces = goodFullFitYint[[All, "PhaseTrace"]];
fullFitYintReTraces    = goodFullFitYint[[All, "ReTrace"]];
fullFitYintImTraces    = goodFullFitYint[[All, "ImTrace"]];
fitYintColors          = ColorData["Rainbow"] /@ Rescale[Range[Length[goodFullFitYint]]];

(* ---------- Residual traces ---------- *)
residualYintResults =
	Table[
		Module[
			{
				rawAssoc, fitAssoc, fRaw, yRaw, yFit,
				magResTrace, phaseResTrace, reResTrace, imResTrace
			},

			rawAssoc = rawYintResults[[k]];
			fitAssoc = fullFitYintResults[[k]];

			If[!TrueQ[rawAssoc["OK"]] || !TrueQ[fitAssoc["OK"]],
				Return[
					<|
						"OK" -> False,
						"Message" -> "raw/fit pair unavailable"
					|>
				]
			];

			fRaw = rawAssoc["FreqHz"];
			yRaw = rawAssoc["YIntRaw"];
			yFit = fitAssoc["YIntFit"];

			If[Length[fRaw] =!= Length[yRaw] || Length[yRaw] =!= Length[yFit],
				Return[
					<|
						"OK" -> False,
						"Message" -> "raw and fit lengths mismatch"
					|>
				]
			];

			magResTrace   = Transpose[{fRaw, Abs[(yRaw - yFit)/yRaw]}];
			phaseResTrace = Transpose[{fRaw, Abs[(180./Pi) Arg[yRaw] - (180./Pi) Arg[yFit]]}];
			reResTrace    = Transpose[{fRaw, Abs[Re[yRaw] - Re[yFit]]}];
			imResTrace    = Transpose[{fRaw, Abs[(-Im[yRaw]) - (-Im[yFit])]}];

			<|
				"OK" -> True,
				"MagResTrace" -> magResTrace,
				"PhaseResTrace" -> phaseResTrace,
				"ReResTrace" -> reResTrace,
				"ImResTrace" -> imResTrace
			|>
		],
		{k, Length[files]}
	];

goodResidualYint = Select[residualYintResults, TrueQ[#["OK"]] &];

If[Length[goodResidualYint] == 0,
	Print["No residual admittance traces could be computed."];
	Abort[];
];

residualYintMagTraces   = goodResidualYint[[All, "MagResTrace"]];
residualYintPhaseTraces = goodResidualYint[[All, "PhaseResTrace"]];
residualYintReTraces    = goodResidualYint[[All, "ReResTrace"]];
residualYintImTraces    = goodResidualYint[[All, "ImResTrace"]];
resYintColors           = ColorData["Rainbow"] /@ Rescale[Range[Length[goodResidualYint]]];

(* ---------- Shared y-axis ranges ---------- *)

allMagY = Select[
	Join[
		Flatten[rawYintMagTraces[[All, All, 2]]],
		Flatten[fullFitYintMagTraces[[All, All, 2]]]
	],
	NumericQ[#] && # > 0 &
];

allPhaseY = Select[
	Join[
		Flatten[rawYintPhaseTraces[[All, All, 2]]],
		Flatten[fullFitYintPhaseTraces[[All, All, 2]]]
	],
	NumericQ
];

allReY = Select[
	Join[
		Flatten[rawYintReTraces[[All, All, 2]]],
		Flatten[fullFitYintReTraces[[All, All, 2]]]
	],
	NumericQ
];

allImY = Select[
	Join[
		Flatten[rawYintImTraces[[All, All, 2]]],
		Flatten[fullFitYintImTraces[[All, All, 2]]]
	],
	NumericQ
];

allMagResY = Select[Flatten[residualYintMagTraces[[All, All, 2]]], NumericQ[#] && # > 0 &];
allPhaseResY = Select[Flatten[residualYintPhaseTraces[[All, All, 2]]], NumericQ];
allReResY = Select[Flatten[residualYintReTraces[[All, All, 2]]], NumericQ];
allImResY = Select[Flatten[residualYintImTraces[[All, All, 2]]], NumericQ];

If[
	allMagY === {} || allPhaseY === {} || allReY === {} || allImY === {} ||
	allMagResY === {} || allPhaseResY === {} || allReResY === {} || allImResY === {},
	Print["Could not determine shared axis ranges."];
	Abort[];
];

magYRange      = {Min[allMagY], Max[allMagY]};
phaseYRange    = {Min[allPhaseY], Max[allPhaseY]};
reYRange       = {Min[allReY], Max[allReY]};
imYRange       = {Min[allImY], Max[allImY]};
magResYRange   = {Min[allMagResY], Max[allMagResY]};
phaseResYRange = {Min[allPhaseResY], Max[allPhaseResY]};
reResYRange    = {Min[allReResY], Max[allReResY]};
imResYRange    = {Min[allImResY], Max[allImResY]};

Print["Shared |Y| y-range = ", magYRange];
Print["Shared phase y-range = ", phaseYRange];
Print["Shared Re(Y) y-range = ", reYRange];
Print["Shared -Im(Y) y-range = ", imYRange];

(* ---------- Plots ---------- *)

rawYintMagPlot =
	ListLogLogPlot[
		rawYintMagTraces,
		PlotStyle -> rawYintColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "|Yint,raw| (S)"},
		PlotRange -> {All, magYRange},
		ImageSize -> 900,
		Background -> White,
		PlotLabel -> Row[{"RAW |Yint|   n=", Length[goodRawYint]}]
	];

fitYintMagPlot =
	ListLogLogPlot[
		fullFitYintMagTraces,
		PlotStyle -> fitYintColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "|Yint,fit| (S)"},
		PlotRange -> {All, magYRange},
		ImageSize -> 900,
		Background -> White,
		PlotLabel -> Row[{"FIT |Yint|   n=", Length[goodFullFitYint]}]
	];

resYintMagPlot =
	ListLogLogPlot[
		residualYintMagTraces,
		PlotStyle -> resYintColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "|Yraw - Yfit| (S)"},
		PlotRange -> {All, magResYRange},
		ImageSize -> 900,
		Background -> White,
		PlotLabel -> "ABS residual |Yraw - Yfit|"
	];

rawYintPhasePlot =
	ListLogLinearPlot[
		rawYintPhaseTraces,
		PlotStyle -> rawYintColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "Phase(Yint,raw) (deg)"},
		PlotRange -> {All, phaseYRange},
		ImageSize -> 900,
		Background -> White,
		PlotLabel -> Row[{"RAW phase(Yint)   n=", Length[goodRawYint]}]
	];

fitYintPhasePlot =
	ListLogLinearPlot[
		fullFitYintPhaseTraces,
		PlotStyle -> fitYintColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "Phase(Yint,fit) (deg)"},
		PlotRange -> {All, phaseYRange},
		ImageSize -> 900,
		Background -> White,
		PlotLabel -> Row[{"FIT phase(Yint)   n=", Length[goodFullFitYint]}]
	];

resYintPhasePlot =
	ListLogLinearPlot[
		residualYintPhaseTraces,
		PlotStyle -> resYintColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "|Phase raw - fit| (deg)"},
		PlotRange -> {All, phaseResYRange},
		ImageSize -> 900,
		Background -> White,
		PlotLabel -> "ABS phase residual"
	];

rawYintRePlot =
	ListLogLinearPlot[
		rawYintReTraces,
		PlotStyle -> rawYintColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "Re(Yint,raw) (S)"},
		PlotRange -> {All, reYRange},
		ImageSize -> 900,
		Background -> White,
		PlotLabel -> "RAW Re(Yint)"
	];

fitYintRePlot =
	ListLogLinearPlot[
		fullFitYintReTraces,
		PlotStyle -> fitYintColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "Re(Yint,fit) (S)"},
		PlotRange -> {All, reYRange},
		ImageSize -> 900,
		Background -> White,
		PlotLabel -> "FIT Re(Yint)"
	];

resYintRePlot =
	ListLogLinearPlot[
		residualYintReTraces,
		PlotStyle -> resYintColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "|Re raw - fit| (S)"},
		(* PlotRange -> {All, reResYRange}, *)
		PlotRange -> {All, {0,0.0001}},
		ImageSize -> 900,
		Background -> White,
		PlotLabel -> "ABS Re residual"
	];

rawYintImPlot =
	ListLogLinearPlot[
		rawYintImTraces,
		PlotStyle -> rawYintColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "-Im(Yint,raw) (S)"},
		PlotRange -> {All, imYRange},
		ImageSize -> 900,
		Background -> White,
		PlotLabel -> "RAW -Im(Yint)"
	];

fitYintImPlot =
	ListLogLinearPlot[
		fullFitYintImTraces,
		PlotStyle -> fitYintColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "-Im(Yint,fit) (S)"},
		PlotRange -> {All, imYRange},
		ImageSize -> 900,
		Background -> White,
		PlotLabel -> "FIT -Im(Yint)"
	];

resYintImPlot =
	ListLogLinearPlot[
		residualYintImTraces,
		PlotStyle -> resYintColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "|(-Im raw) - (-Im fit)| (S)"},
		(* PlotRange -> {All, imResYRange}, *)
		PlotRange -> {All, {0,0.0001}},
		ImageSize -> 900,
		Background -> White,
		PlotLabel -> "ABS -Im residual"
	];

GraphicsGrid[
	{
		{rawYintMagPlot, fitYintMagPlot, resYintMagPlot},
		{rawYintPhasePlot, fitYintPhasePlot, resYintPhasePlot},
		{rawYintRePlot, fitYintRePlot, resYintRePlot},
		{rawYintImPlot, fitYintImPlot, resYintImPlot}
	}, ImageSize->1200
]


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ============================== *)
(* User control: area window      *)
(* ============================== *)
fMinArea = 30/(2 Pi);
fMaxArea = 500/(2 Pi);
fMid = 150/(2 Pi);   (* Hz, user-selected split point *)

kMinArea = 2 Pi fMinArea;   (* s^-1 *)
kMaxArea = 2 Pi fMaxArea;   (* s^-1 *)
kMid = 2 Pi fMid;           (* s^-1 *)

(* ============================== *)
(* Helper: integrate g over log10(k) *)
(* g is in F/decade, so integrate vs log10(k) *)
(* ============================== *)

areaOverLog10k[k_List, g_List] := Module[
  {ord, kk, gg, x},
  If[Length[k] < 2 || Length[g] < 2, Return[0.0]];
  ord = Ordering[k];
  kk = k[[ord]];
  gg = g[[ord]];
  x = Log10[kk];
  N @ Total[Differences[x] * MovingAverage[gg, 2]]
];

(* ============================== *)
(* Compute areas for each dataset *)
(* restricted to kMinArea <= k <= kMaxArea *)
(* then split at kMid             *)
(* ============================== *)

areaResults = Table[
  Module[
    {
      spec, k, g, ord, kk, gg, keep,
      kkUse, ggUse, idxSplit, cSlow, cFast, cTot
    },

    spec = goodSpecs[[i, "Spec"]];
    k = 1/spec["Tau"];
    g = spec["g"];

    (* sort by k increasing *)
    ord = Ordering[k];
    kk = k[[ord]];
    gg = g[[ord]];

    (* restrict to requested area window *)
    keep = (kMinArea <= # <= kMaxArea) & /@ kk;
    kkUse = Pick[kk, keep];
    ggUse = Pick[gg, keep];

    (* split inside restricted window *)
    idxSplit = LengthWhile[kkUse, # <= kMid &];

    cSlow =
      If[idxSplit >= 2,
        areaOverLog10k[kkUse[[;; idxSplit]], ggUse[[;; idxSplit]]],
        0.0
      ];

    cFast =
      If[idxSplit + 1 <= Length[kkUse] - 1,
        areaOverLog10k[kkUse[[idxSplit + 1 ;;]], ggUse[[idxSplit + 1 ;;]]],
        0.0
      ];

    cTot = cSlow + cFast;

    <|
      "i" -> i,
      "File" -> goodSpecs[[i, "File"]],
      "fMinArea_Hz" -> fMinArea,
      "fMid_Hz" -> fMid,
      "fMaxArea_Hz" -> fMaxArea,
      "kMinArea_s^-1" -> kMinArea,
      "kMid_s^-1" -> kMid,
      "kMaxArea_s^-1" -> kMaxArea,
      "cSlow_F" -> cSlow,   (* kMinArea <= k <= kMid *)
      "cFast_F" -> cFast,   (* kMid < k <= kMaxArea *)
      "cTotal_F" -> cTot,
      "FractionBound" -> If[cTot > 0, cFast/cTot, Indeterminate]
    |>
  ],
  {i, Length[goodSpecs]}
];

(* Quick look *)
Dataset[areaResults];
(* ============================== *)
(* User control: Langmuir KD       *)
(* ============================== *)
KD = 250;
KD = 144;  (* user-defined; same concentration units you want out, e.g. uM *)
(* 144 uM was the published value measured with EIS https://pubs.acs.org/doi/full/10.1021/acssensors.3c00632 *)

(* ============================== *)
(* Helper: Langmuir inversion      *)
(* ============================== *)
langmuirConcFromF[f_?NumericQ, kd_?NumericQ] := Module[{eps = 10^-12, ff},
  (* clip f away from exactly 0 or 1 to avoid division blowups *)
  ff = Clip[f, {0 + eps, 1 - eps}];
  kd * ff/(1 - ff)
];

(* ============================== *)
(* Add concentration to each row   *)
(* ============================== *)

areaResultsWithConc =
  Map[
    Function[assoc,
      Module[{f = assoc["FractionBound"], cEst},
        cEst = langmuirConcFromF[f, KD];
        Join[assoc, <|
          "KD" -> KD,
          "Conc_Est" -> cEst
        |>]
      ]
    ],
    areaResults
  ];

Dataset[areaResultsWithConc];

(* ================================= *)
(* Extract number before "uM"        *)
(* ================================= *)

getConcFromFile[file_String] := Module[
  {name, hit},
  
  name = FileNameTake[file];
  
  hit = StringCases[
    name,
    RegularExpression["(\\d+(?:\\.\\d+)?)uM"] :> "$1"
  ];
  
  If[hit === {},
    Missing["NoMatch"],
    ToExpression[First[hit]]
  ]
];

(* ================================= *)
(* Add column to areaResults         *)
(* ================================= *)

areaResultsWithConc =
  Map[
    Function[assoc,
      Join[
        assoc,
        <|"Conc_Actual" -> getConcFromFile[assoc["File"]]|>
      ]
    ],
    areaResults
  ];

Dataset[areaResultsWithConc];
(* ============================== *)
(* Combine actual + estimated conc *)
(* ============================== *)

areaResultsWithConc =
  Map[
    Function[assoc,
      Module[{f, cEst, cAct},
        f = assoc["FractionBound"];
        cEst = langmuirConcFromF[f, KD];
        cAct = getConcFromFile[assoc["File"]];

        Join[assoc, <|
          "KD" -> KD,
          "Conc_Actual" -> cAct,
          "Conc_Est" -> cEst
        |>]
      ]
    ],
    areaResults
  ];

Dataset[areaResultsWithConc];



concTable = areaResultsWithConc[[All, {"File", "Conc_Actual", "Conc_Est"}]];
Dataset[concTable];
(* ============================================ *)
(* 1) Helpers: parse concentration + electrode  *)
(* ============================================ *)

getConcFromFile[file_String] := Module[{name, hit},
  name = FileNameTake[file];
  hit = StringCases[
    name,
    RegularExpression["(?i)(\\d+(?:\\.\\d+)?)\\s*uM"] :> "$1"
  ];
  If[hit === {}, Missing["NoConc"], ToExpression[First[hit]]]
];

getElectrodeFromFile[file_String] := Module[{name, hit},
  name = FileNameTake[file];
  hit = StringCases[
    name,
    RegularExpression["E(\\d+)"] :> ("E" <> "$1")
  ];
  If[hit === {}, Missing["NoElectrode"], First[hit]]
];

(* ============================================ *)
(* 2) Langmuir inversion                        *)
(* ============================================ *)
langmuirConcFromF[f_?NumericQ, kd_?NumericQ] := Module[{eps = 10^-12, ff},
  ff = Clip[f, {eps, 1 - eps}];
  kd * ff/(1 - ff)
];

(* ============================================ *)
(* 3) Build one combined table                  *)
(* ============================================ *)

areaResultsFull = Table[
  Module[{spec, file, k, g, ord, kk, gg, keep, kkUse, ggUse, idx, cSlow, cFast, cTot, frac},
    spec = goodSpecs[[i, "Spec"]];
    file = goodSpecs[[i, "File"]];
    k = 1/spec["Tau"];
    g = spec["g"];

    ord = Ordering[k];
    kk = k[[ord]];
    gg = g[[ord]];

    keep = (kMinArea <= # <= kMaxArea) & /@ kk;
    kkUse = Pick[kk, keep];
    ggUse = Pick[gg, keep];

    idx = LengthWhile[kkUse, # <= kMid &];

    cSlow = If[idx >= 2, areaOverLog10k[kkUse[[;; idx]], ggUse[[;; idx]]], 0.0];
    cFast = If[idx + 1 <= Length[kkUse] - 1, areaOverLog10k[kkUse[[idx + 1 ;;]], ggUse[[idx + 1 ;;]]], 0.0];
    cTot = cSlow + cFast;
    frac = If[cTot > 0, cFast/cTot, Indeterminate];

    <|
      "i" -> i,
      "TimeHr" -> expTimesHr[[i]],
      "File" -> file,
      "Electrode" -> getElectrodeFromFile[file],
      "ConcActual" -> getConcFromFile[file],
      "KD" -> KD,
      "cSlow_F" -> cSlow,
      "cFast_F" -> cFast,
      "cTotal_F" -> cTot,
      "FractionBound" -> frac,
      "ConcEst" -> If[NumericQ[frac], langmuirConcFromF[frac, KD], Missing["NoFrac"]]
    |>
  ],
  {i, Length[goodSpecs]}
];
Dataset[areaResultsFull];
(* ============================================ *)
(* Build plotting table from areaResultsWithConc *)
(* Requires keys: "Conc_Actual", "Conc_Est", "File" *)
(* ============================================ *)
(* ============================================ *)
(* Build plotting table                          *)
(* areaResultsFull must contain: "Conc_Actual", "Conc_Est", "File" *)
(* ============================================ *)

plotRows =
  Select[
    Map[
      Function[assoc,
        Module[{x, y, e},
          x = assoc["ConcActual"];
          y = assoc["ConcEst"];
          e = getElectrodeFromFile[assoc["File"]];
          <|
            "ConcActual" -> x,
            "ConcEst" -> y,
            "Electrode" -> e,
            "File" -> assoc["File"]
          |>
        ]
      ],
      areaResultsFull
    ],
    (NumericQ[#["ConcActual"]] && NumericQ[#["ConcEst"]] && StringQ[#["Electrode"]]) &
  ];

Dataset[plotRows];
(* ============================================ *)
(* Plot: Conc_Est vs Conc_Actual, colored by electrode *)
(* ============================================ *)
groups = GroupBy[plotRows, #["Electrode"] &];

traces =
  AssociationMap[
    ({#["ConcActual"], #["ConcEst"]} & /@ groups[#]) &,
    Keys[groups]
  ];

traces;

colors = <|
  "E1" -> Red,
  "E2" -> Blue,
  "E3" -> Darker[Green]
|>;

(* ============================================ *)
(* Build plotting table: FractionBound vs Actual *)
(* Requires: areaResultsFull has keys:
     "Conc_Actual", "FractionBound", "File"
   and getElectrodeFromFile[file] returns "E1"/"E2"/"E3"
*)
(* ============================================ *)

plotRowsFB =
  Select[
    Map[
      Function[assoc,
        Module[{x, f, e},
          x = assoc["ConcActual"];
          f = assoc["FractionBound"];
          e = getElectrodeFromFile[assoc["File"]];
          <|
            "ConcActual" -> x,
            "FractionBound" -> f,
            "Electrode" -> e,
            "File" -> assoc["File"]
          |>
        ]
      ],
      areaResultsFull
    ],
    (NumericQ[#["ConcActual"]] &&
     NumericQ[#["FractionBound"]] &&
     StringQ[#["Electrode"]]) &
  ];

Dataset[plotRowsFB];
tracesFB = GroupBy[
  plotRowsFB,
  #["Electrode"] &,
  ( {#["ConcActual"], #["FractionBound"]} & /@ # ) &
];

(* ============================================ *)
(* Plot: FractionBound vs Actual Concentration  *)
(* + overlay Langmuir prediction                *)
(* ============================================ *)

colors = <|
  "E1" -> Red,
  "E2" -> Blue,
  "E3" -> Darker[Green]
|>;

(* Domain for the Langmuir curve: based on your data *)
allConc = Flatten[Values[tracesFB], 1][[All, 1]];
cMin = Max[100, Min[Select[allConc, NumericQ]]];   (* avoid 0 on log scale *)
cMax = Max[Select[allConc, NumericQ]];

langmuirFB[c_?NumericQ, kd_?NumericQ] := c/(c + kd);

Show[
  {
    ListPlot[
      Values[tracesFB],
      PlotStyle -> (colors /@ Keys[tracesFB]),
      PlotMarkers -> {Automatic, 9},
      Joined -> False,
      PlotRange -> {{0, 5 cMax}, {0, 1}}
    ],
    Plot[
      langmuirFB[c, KD],
      {c, 0, 5 cMax},
      PlotStyle -> {GrayLevel[0.35], Thick}
    ],
    Plot[
      langmuirFB[c, KD] + 0.08,
      {c, 0, 5 cMax},
      PlotStyle -> {GrayLevel[0.35], Thick}
    ]
  },
  Frame -> True,
  FrameLabel -> {"Actual concentration (uM)", "Fraction bound"},
  PlotLegends -> Placed[
    LineLegend[
      Join[
        (Style[#, colors[#]] & /@ Keys[tracesFB]),
        {Style["Langmuir (KD=" <> ToString[KD] <> " uM)", GrayLevel[0.35]]}
      ]
    ],
    Right
  ],
    ImageSize -> 700,
  Background -> White
]



(* ---------- Fraction bound vs time ---------- *)

timeSeriesRows = SortBy[
  Select[
    areaResultsFull,
    NumericQ[#["FractionBound"]] && IntegerQ[#["i"]] &
  ],
  #["i"] &
];

fbTimeData = Transpose[{timeSeriesRows[[All, "i"]], timeSeriesRows[[All, "FractionBound"]]}];

fbTimeRows = SortBy[
  Select[
    areaResultsFull,
    NumericQ[#["TimeHr"]] && NumericQ[#["FractionBound"]] &
  ],
  #["TimeHr"] &
];

fbTimeData = ({#["TimeHr"], #["FractionBound"]} & /@ fbTimeRows);

ListPlot[
  fbTimeData,
  Joined -> True,
  PlotMarkers -> {Automatic, 8},
  Background -> White,
  PlotRangePadding -> Scaled[0.02],
  Frame -> True,
  FrameStyle -> Directive[Black, AbsoluteThickness[1.2]],
  FrameLabel -> {
    Style["Time (hours)", 18, Black],
    Style["Fraction folded", 18, Black]
  },
  BaseStyle -> {FontFamily -> "Arial", 14, Black},
  LabelStyle -> Directive[18, Black],
  FrameTicksStyle -> Directive[14, Black],
  FrameTicks -> {
    {Automatic, None},
    {Automatic, None}
  },
  PlotRange -> {All, {0, 1}},
  ImageSize -> 900,
  PlotLabel -> Style["Fraction folded vs time", 16, Black]
]


(* ============================================================ *)
(* TIME-SERIES / PSD ANALYSIS OF EIS FOURIER COMPONENTS         *)
(* Uses already-built goodRawZ and goodSpecs                    *)
(* Last 10 hours only                                           *)
(* Outputs:                                                     *)
(*   1) Overlay PSD plot for all EIS frequencies               *)
(*   2) Table of mean Z quantities vs frequency                *)
(*   3) Variance vs mean for impedance                         *)
(*   4) Variance vs mean for admittance after subtracting Rs   *)
(* ============================================================ *)

If[Length[goodRawZ] == 0 || Length[goodSpecs] == 0,
	Print["goodRawZ or goodSpecs is empty."];
	Abort[];
];

(* ---------- helper to pull finish times ---------- *)
finishTimesS = Lookup[goodSpecs[[All, "Spec"]], "FinishTimeS", Missing["NoTime"]];

validTimeIdx = Select[
	Range[Length[goodRawZ]],
	NumericQ[finishTimesS[[#]]] &
];

If[Length[validTimeIdx] < 20,
	Print["Too few spectra with valid FinishTimeS."];
	Abort[];
];

timesSAll = finishTimesS[[validTimeIdx]];
rawSubsetAll = goodRawZ[[validTimeIdx]];
specSubsetAll = goodSpecs[[validTimeIdx]];

tMax = Max[timesSAll];
tMinKeep = tMax - 20.*3600.;

keepLast10h = Map[# >= tMinKeep &, timesSAll];

timesS = Pick[timesSAll, keepLast10h];
rawSubset = Pick[rawSubsetAll, keepLast10h];
specSubset = Pick[specSubsetAll, keepLast10h];

If[Length[rawSubset] < 20,
	Print["Too few spectra in last 20 hours."];
	Abort[];
];

(* ---------- sort by time ---------- *)
ord = Ordering[timesS];
timesS = timesS[[ord]];
rawSubset = rawSubset[[ord]];
specSubset = specSubset[[ord]];

timesHr = (timesS - Min[timesS])/3600.0;
nSpec = Length[rawSubset];

Print["Using ", nSpec, " spectra from last 10 hours."];

(* ---------- common frequency grid check ---------- *)
freqTemplate = rawSubset[[1, "FreqHz"]];

sameFreqGridQ = And @@ Table[
	rawSubset[[k, "FreqHz"]] === freqTemplate,
	{k, 2, nSpec}
];

If[!sameFreqGridQ,
	Print["Frequency grids are not identical across spectra."];
	Abort[];
];

freqList = N[freqTemplate];
nFreq = Length[freqList];

(* ---------- build matrices ---------- *)
zMat = Table[
	rawSubset[[k, "ZRaw"]],
	{k, 1, nSpec}
];

rsVec = Table[
	specSubset[[k, "Spec", "Rs"]],
	{k, 1, nSpec}
];

If[!VectorQ[rsVec, NumericQ],
	Print["Could not extract numeric Rs values from goodSpecs."];
	Abort[];
];

yMat = Table[
	1/(zMat[[k]] - rsVec[[k]]),
	{k, 1, nSpec}
];

(* ---------- sampling interval estimate ---------- *)
dtList = Differences[timesS];
dtMed = Median[dtList];

If[!NumericQ[dtMed] || dtMed <= 0,
	Print["Could not determine positive sampling interval."];
	Abort[];
];

fsTime = 1./dtMed;
Print["Median sampling interval = ", NumberForm[dtMed, {8, 2}], " s"];
Print["Effective time-series sample rate = ", NumberForm[fsTime, {8, 5}], " Hz"];

(* ---------- simple periodogram PSD helper ---------- *)
clearPSD[x_List, dt_?NumericQ] := Module[
	{
		x0, n, fs, fft, kmax, psd, ff
	},
	n = Length[x];
	If[n < 4, Return[<|"PSDfreq" -> {}, "PSD" -> {}|>]];

	fs = 1./dt;
	x0 = N[x - Mean[x]];
	fft = Fourier[x0, FourierParameters -> {1, -1}];

	kmax = Floor[n/2];

	ff = Range[0, kmax] * fs/n;
	psd = (1./(fs*n)) * Abs[fft[[1 ;; kmax + 1]]]^2;

	If[kmax >= 1,
		psd[[2 ;; -2]] = 2.0 * psd[[2 ;; -2]];
	];

	<|"PSDfreq" -> ff, "PSD" -> psd|>
];

(* ---------- compute PSD and statistics at each EIS frequency ---------- *)
psdResults = Table[
	Module[
		{
			f = freqList[[j]],
			zSeries, ySeries,
			zMagSeries, zReSeries, zImSeries,
			yMagSeries, yReSeries, yImSeries,
			psdZMag, meanZMag, varZMagRaw, varZMagPSD,
			psdYMag, meanYMag, varYMagRaw, varYMagPSD,
			meanZRe, meanNegZIm, varZReRaw, varNegZImRaw,
			meanYRe, meanNegYIm, varYReRaw, varNegYImRaw
		},

		zSeries = zMat[[All, j]];
		ySeries = yMat[[All, j]];

		zMagSeries = Abs[zSeries];
		zReSeries = Re[zSeries];
		zImSeries = -Im[zSeries];

		yMagSeries = Abs[ySeries];
		yReSeries = Re[ySeries];
		yImSeries = -Im[ySeries];

		psdZMag = clearPSD[zMagSeries, dtMed];
		psdYMag = clearPSD[yMagSeries, dtMed];

		meanZMag = Mean[zMagSeries];
		varZMagRaw = Variance[zMagSeries];
		varZMagPSD = Total[Most[psdZMag["PSD"]] * Differences[psdZMag["PSDfreq"]]];

		meanYMag = Mean[yMagSeries];
		varYMagRaw = Variance[yMagSeries];
		varYMagPSD = Total[Most[psdYMag["PSD"]] * Differences[psdYMag["PSDfreq"]]];

		meanZRe = Mean[zReSeries];
		meanNegZIm = Mean[zImSeries];
		varZReRaw = Variance[zReSeries];
		varNegZImRaw = Variance[zImSeries];

		meanYRe = Mean[yReSeries];
		meanNegYIm = Mean[yImSeries];
		varYReRaw = Variance[yReSeries];
		varNegYImRaw = Variance[yImSeries];

		<|
			"Freq" -> f,

			"PSDfreq" -> psdZMag["PSDfreq"],
			"PSD" -> psdZMag["PSD"],

			"MeanZMag" -> meanZMag,
			"VarZMagRaw" -> varZMagRaw,
			"VarZMagPSD" -> varZMagPSD,
			"MeanZRe" -> meanZRe,
			"MeanNegZIm" -> meanNegZIm,
			"VarZReRaw" -> varZReRaw,
			"VarNegZImRaw" -> varNegZImRaw,

			"MeanYMag" -> meanYMag,
			"VarYMagRaw" -> varYMagRaw,
			"VarYMagPSD" -> varYMagPSD,
			"MeanYRe" -> meanYRe,
			"MeanNegYIm" -> meanNegYIm,
			"VarYReRaw" -> varYReRaw,
			"VarNegYImRaw" -> varNegYImRaw
		|>
	],
	{j, 1, nFreq}
];

goodPSDResults = Select[
	psdResults,
	AssociationQ[#] &&
	KeyExistsQ[#, "PSDfreq"] &&
	KeyExistsQ[#, "PSD"] &&
	VectorQ[#["PSDfreq"], NumericQ] &&
	VectorQ[#["PSD"], NumericQ] &&
	Length[#["PSDfreq"]] == Length[#["PSD"]] &&
	Length[#["PSDfreq"]] > 2 &
];

If[Length[goodPSDResults] == 0,
	Print["No valid PSD results found."];
	Abort[];
];

(* ---------- overlay PSD plot ---------- *)
psdCurves = Transpose[{#["PSDfreq"]*3600, #["PSD"]}] & /@ goodPSDResults;

psdLabels = Table[
	ToString@NumberForm[goodPSDResults[[k, "Freq"]], {7, 2}] <> " Hz",
	{k, Length[goodPSDResults]}
];

psdColors = ColorData["Rainbow"] /@ Rescale[Range[Length[psdCurves]]];

psdOverlayPlot =
	ListLogLogPlot[
		psdCurves,
		PlotStyle -> psdColors,
		Frame -> True,
		Joined->True,
		Axes -> False,
		FrameLabel -> {"Temporal frequency (Hz)", "PSD of |Z|"},
		PlotRange -> All,
		ImageSize -> 900,
		Background -> White,
		PlotLegends -> Placed[
			LineLegend[
				psdColors,
				psdLabels,
				LegendLayout -> "Column",
				LabelStyle -> Directive[11]
			],
			Right
		],
		PlotLabel -> "Overlay PSD of |Z| for each EIS Fourier component"
	];
psdOverlayPlot




(* ---------- table of means ---------- *)
meanTable = Table[
	{
		psdResults[[j, "Freq"]],
		psdResults[[j, "MeanZMag"]],
		psdResults[[j, "MeanZRe"]],
		psdResults[[j, "MeanNegZIm"]],
		psdResults[[j, "MeanYMag"]],
		psdResults[[j, "MeanYRe"]],
		psdResults[[j, "MeanNegYIm"]]
	},
	{j, 1, Length[psdResults]}
];

meanTableHeader = {
	"Freq (Hz)",
	"Mean |Z|",
	"Mean Re[Z]",
	"Mean -Im[Z]",
	"Mean |Y|",
	"Mean Re[Y]",
	"Mean -Im[Y]"
};

meanTableGrid =
	Grid[
		Prepend[
			NumberForm[#, {10, 4}] & /@ meanTable,
			meanTableHeader
		],
		Frame -> All,
		Background -> {None, {White}},
		ItemStyle -> Directive[12]
	];

(* ---------- variance vs mean traces ---------- *)
zVarMeanRaw = Select[
	Table[
		{psdResults[[j, "MeanZMag"]], psdResults[[j, "VarZMagRaw"]]},
		{j, 1, Length[psdResults]}
	],
	NumericQ[#[[1]]] && NumericQ[#[[2]]] && #[[1]] > 0 && #[[2]] > 0 &
];

zVarMeanPSD = Select[
	Table[
		{psdResults[[j, "MeanZMag"]], psdResults[[j, "VarZMagPSD"]]},
		{j, 1, Length[psdResults]}
	],
	NumericQ[#[[1]]] && NumericQ[#[[2]]] && #[[1]] > 0 && #[[2]] > 0 &
];

yVarMeanRaw = Select[
	Table[
		{psdResults[[j, "MeanYMag"]], psdResults[[j, "VarYMagRaw"]]},
		{j, 1, Length[psdResults]}
	],
	NumericQ[#[[1]]] && NumericQ[#[[2]]] && #[[1]] > 0 && #[[2]] > 0 &
];

yVarMeanPSD = Select[
	Table[
		{psdResults[[j, "MeanYMag"]], psdResults[[j, "VarYMagPSD"]]},
		{j, 1, Length[psdResults]}
	],
	NumericQ[#[[1]]] && NumericQ[#[[2]]] && #[[1]] > 0 && #[[2]] > 0 &
];

(* ---------- reference scaling lines ---------- *)
makeReferenceLine[data_List, slope_?NumericQ] := Module[
	{x0, y0, a, xs},
	If[Length[data] < 2, Return[{}]];
	x0 = data[[Ceiling[Length[data]/2], 1]];
(* 	y0 = data[[Ceiling[Length[data]/2], 2]]; *);
	y0 = 0.4 data[[Ceiling[Length[data]/2], 2]];

	a = y0/(x0^slope);
	xs = {Min[data[[All, 1]]], Max[data[[All, 1]]]};
	{{xs[[1]], a*xs[[1]]^slope}, {xs[[2]], a*xs[[2]]^slope}}
];

zShotLine = makeReferenceLine[zVarMeanRaw, 3];

yShotLine = makeReferenceLine[yVarMeanRaw, 1];

(* ---------- label helper ---------- *)
targetFreqs = {1., 10., 100., 1000., 10000.};

closestIndexToFreq[f_] := First @ Ordering[Abs[N[psdResults[[All, "Freq"]]] - f], 1];

zLabelData = Table[
	Module[{idx, f, x, y},
		idx = closestIndexToFreq[targetFreqs[[k]]];
		f = psdResults[[idx, "Freq"]];
		x = psdResults[[idx, "MeanZMag"]];
		y = psdResults[[idx, "VarZMagRaw"]];
		{f, x, y}
	],
	{k, 1, Length[targetFreqs]}
];

yLabelData = Table[
	Module[{idx, f, x, y},
		idx = closestIndexToFreq[targetFreqs[[k]]];
		f = psdResults[[idx, "Freq"]];
		x = psdResults[[idx, "MeanYMag"]];
		y = psdResults[[idx, "VarYMagRaw"]];
		{f, x, y}
	],
	{k, 1, Length[targetFreqs]}
];

zLabelData = Select[
	zLabelData,
	NumericQ[#[[2]]] && NumericQ[#[[3]]] && #[[2]] > 0 && #[[3]] > 0 &
];

yLabelData = Select[
	yLabelData,
	NumericQ[#[[2]]] && NumericQ[#[[3]]] && #[[2]] > 0 && #[[3]] > 0 &
];

zLabelPrimitives = Flatten@Table[
	With[
		{
			f = zLabelData[[k, 1]],
			pt = zLabelData[[k, {2, 3}]],
			offset = {{1.2, 1.0}, {1.2, -1.0}, {-1.2, 1.0}, {-1.2, -1.0}, {1.2, 0.0}}[[k]]
		},
		{
			Black,
			PointSize[0.018],
			Point[pt],
			Text[
				Style[ToString@NumberForm[f, {6, 0}] <> " Hz", 12, Black, Background -> White],
				pt,
				offset
			]
		}
	],
	{k, 1, Length[zLabelData]}
];

yLabelPrimitives = Flatten@Table[
	With[
		{
			f = yLabelData[[k, 1]],
			pt = yLabelData[[k, {2, 3}]],
			offset = {{1.2, 1.0}, {1.2, -1.0}, {-1.2, 1.0}, {-1.2, -1.0}, {1.2, 0.0}}[[k]]
		},
		{
			Black,
			PointSize[0.018],
			Point[pt],
			Text[
				Style[ToString@NumberForm[f, {6, 0}] <> " Hz", 12, Black, Background -> White],
				pt,
				offset
			]
		}
	],
	{k, 1, Length[yLabelData]}
];

zVarPlot =
	ListLogLogPlot[
		{
			zVarMeanRaw,
			zVarMeanPSD,
			zShotLine,
			zAltLine
		},
		Joined -> {False, False, True, True},
		PlotMarkers -> {
			{Automatic, Medium},
			{Automatic, Medium},
			None,
			None
		},
		PlotStyle -> {
			Directive[Black],
			Directive[Red],
			Directive[Blue, Dashed, Thick],
			Directive[Darker[Green], Dashed, Thick]
		},
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Mean |Z|", "Variance of |Z|"},
		PlotRange -> All,
		ImageSize -> 800,
		Background -> White,
		PlotLegends -> Placed[
			{
				"Raw variance",
				"PSD-integrated variance",
				"Shot-noise expectation  Var[Z] \[Proportional] |Z|^3",
				"Alternative  Var[Z] \[Proportional] |Z|^1"
			},
			Right
		],
		PlotLabel -> "Impedance variance vs mean",
		Epilog -> zLabelPrimitives
	];

yVarPlot =
	ListLogLogPlot[
		{
			yVarMeanRaw,
			yVarMeanPSD,
			yShotLine,
			yAltLine
		},
		Joined -> {False, False, True, True},
		PlotMarkers -> {
			{Automatic, Medium},
			{Automatic, Medium},
			None,
			None
		},
		PlotStyle -> {
			Directive[Black],
			Directive[Red],
			Directive[Blue, Dashed, Thick],
			Directive[Darker[Green], Dashed, Thick]
		},
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Mean |Y|", "Variance of |Y|"},
		PlotRange -> All,
		ImageSize -> 800,
		Background -> White,
		PlotLegends -> Placed[
			{
				"Raw variance",
				"PSD-integrated variance",
				"Shot-noise expectation  Var[Y] \[Proportional] |Y|^1",
				"Alternative  Var[Y] \[Proportional] |Y|^2"
			},
			Right
		],
		PlotLabel -> "Admittance variance vs mean",
		Epilog -> yLabelPrimitives
	];

(* ---------- display ---------- *)

meanTableGrid;
zVarPlot;
yVarPlot;


(* ------------------------------------------------------------ *)
(* Build labeled subset of Z variance data                      *)
(* ------------------------------------------------------------ *)

targetFreqs = {1, 10, 100, 300,500,1000, 10000};

(* extract mean-Z vs variance pairs with their frequencies *)
zVarFull = Table[
	{
		psdResults[[j, "Freq"]],
		psdResults[[j, "MeanZMag"]],
		psdResults[[j, "VarZMagRaw"]]
	},
	{j, Length[psdResults]}
];

zVarFull = Select[
	zVarFull,
	NumericQ[#[[2]]] && NumericQ[#[[3]]] && #[[2]] > 0 && #[[3]] > 0 &
];

(* find closest frequency rows *)
labelRows =
	Table[
		First @ MinimalBy[zVarFull, Abs[#[[1]] - f] &],
		{f, targetFreqs}
	];

(* convert to Callout points *)
zLabelPoints =
	Table[
		Callout[
			{row[[2]], row[[3]]},
			Style[ToString[Round[row[[1]]]] <> " Hz", 14, Black],
			Above,
			Appearance -> "Leader"
		],
		{row, labelRows}
	];

(* unlabeled base data *)
zBasePoints = zVarFull[[All, {2, 3}]];

(* ------------------------------------------------------------ *)
(* Final impedance variance plot with labels                    *)
(* ------------------------------------------------------------ *)

zVarPlot =
	ListLogLogPlot[
		{
			zBasePoints,
			zLabelPoints,
			zShotLine
		},
		Joined -> {False, False, True, True},
		PlotMarkers -> {
			{Automatic, Medium},
			None,
			None
		},
		PlotStyle -> {
			Directive[Black],
			Directive[Black],
			Directive[Blue, Dashed, Thick]
		},
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Mean |Z|", "Variance of |Z|"},
		PlotRange -> All,
		ImageSize -> 800,
		Background -> White,
		PlotLegends -> Placed[
			{
				"Raw variance",
				"Labeled representative frequencies",
				"Shot noise slope"
			},
			Right
		],
		PlotLabel -> "Impedance variance vs mean (labeled)"
	];



(* ------------------------------------------------------------ *)
(* Build labeled subset of Y variance data                      *)
(* ------------------------------------------------------------ *)


(* extract mean-Y vs variance pairs with their frequencies *)
yVarFull = Table[
	{
		psdResults[[j, "Freq"]],
		psdResults[[j, "MeanYMag"]],
		psdResults[[j, "VarYMagRaw"]]
	},
	{j, Length[psdResults]}
];

yVarFull = Select[
	yVarFull,
	NumericQ[#[[2]]] && NumericQ[#[[3]]] && #[[2]] > 0 && #[[3]] > 0 &
];

(* find closest frequency rows *)
labelRowsY =
	Table[
		First @ MinimalBy[yVarFull, Abs[#[[1]] - f] &],
		{f, targetFreqs}
	];

(* convert to Callout points *)
yLabelPoints =
	Table[
		Callout[
			{row[[2]], row[[3]]},
			Style[ToString[Round[row[[1]]]] <> " Hz", 14, Black],
			Above,
			Appearance -> "Leader"
		],
		{row, labelRowsY}
	];

(* unlabeled base data *)
yBasePoints = yVarFull[[All, {2, 3}]];

(* ------------------------------------------------------------ *)
(* Final admittance variance plot with labels                   *)
(* ------------------------------------------------------------ *)

yVarPlot =
	ListLogLogPlot[
		{
			yBasePoints,
			yShotLine,
			yLabelPoints
		},
		Joined -> {False,  True,  False},
		PlotMarkers -> {
			{Automatic, Medium},
			None,
			None
		},
		PlotStyle -> {
			Directive[Black],
			Directive[Blue, Dashed, Thick],
			Directive[Black]
		},
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Mean |Y|", "Variance of |Y|"},
		PlotRange -> All,
		ImageSize -> 800,
		Background -> White,
		PlotLegends -> Placed[
			{
				"Raw variance",
				"Shot noise slope",
				"Labeled representative frequencies"
			},
			Right
		],
		PlotLabel -> "Admittance variance vs mean (labeled)"
	];

zVarPlot
yVarPlot





 
