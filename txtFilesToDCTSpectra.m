(* ::Package:: *)

ClearAll["Global`*"]


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


debugFlag = False;
fileDecimation = 1;   (* keep every Nth file: 10 -> ~450/10 = 45 files *)
constantPhaseElementFlag = False;


baseDir = FileNameJoin[{NotebookDirectory[], "data"}] // ExpandFileName


(* examples of switching datasets *)

dataDir = FileNameJoin[{baseDir, "2026-02-25"}];
dataDir = FileNameJoin[{baseDir, "2026-02-25-titration"}]; (* initial titration experiment in PBS *)




dataDir = FileNameJoin[{baseDir, "rebodeplotfcc11sh"}];(* ferrocene data in PBS *)



dataDir = FileNameJoin[{baseDir, "2026-04-13"}]; (* simulated data generated in Mathematica! *)


dataDir = FileNameJoin[{baseDir, "simulated"}];(* simulated data generated in Mathematica! *)
dataDir = FileNameJoin[{baseDir, "2026-04-13"}];(* ferrocene data in PBS take 2 *)



dataDir = FileNameJoin[{baseDir, "spike"}];(* second titration in PBS *)




dataDir = FileNameJoin[{baseDir, "drift_03252026", "E3"}];(* PBS *)


dataDir = FileNameJoin[{baseDir, "04062026_drift\\drift", "E2"}] (* bovine blood *)

dataDir = FileNameJoin[{baseDir, "260507_temp dependence"}] (* ferrocene temperature *)

dataDir = FileNameJoin[{baseDir, "20260514_Fc_drift\\drift", "E1"}] (* ferrocene drift *)



lambdaND = 1 10^-2;

binsPerDecade = 10;
weightPower =-1.0; (* weight each point equally *)

weightPower = -1.25; (* weight each point inverse to variance *)

fMinUse =0.01;      (* Hz *)
fMaxUse = 1000;     (* Hz *)



(* Optional: restrict the frequency range (must match package options) *)
paddingDecades = 0.0; (* sets how many decades beyond the measured frequency range the tau values extend *)
(* 0 = no padding, 1 = a decade of padding. *)

tauMinFactor = 10^-paddingDecades;
tauMaxFactor = 10^paddingDecades;          (* tauMaxUse = tauMax * factor *) (* add factor of 10 to the maximum as padding *)



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
(* 
Print["Files used:"];
Print /@ txtFiles; *)


(* ============================================================ *)
(* 2) FIT EACH FILE TO A DCT SPECTRUM                           *)
(*    - sort by numeric index in parentheses: ... (n).txt       *)
(* ============================================================ *)

fileOrderKey[path_String] := Module[{base, n},
  base = FileBaseName[path];
  n = Quiet @ Check[
     ToExpression @ 
      StringReplace[base, RegularExpression[".*\\((\\d+)\\)$"] -> "$1"],
     Missing["NoIndex"]
   ];
  If[
   IntegerQ[n],
   {StringReplace[base, RegularExpression["\\(\\d+\\)$"] -> ""], n},
   {base, Infinity}
  ]
];

txtFiles = SortBy[txtFiles, fileOrderKey];

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
         "WeightPower" -> weightPower
        ],
        $Failed
       ],
       {NMinimize::dinfeas}
      ];

     ok =
      AssociationQ[spec] &&
       And @@ (KeyExistsQ[spec, #] & /@ {
           "Tau", "g", "C0", "Rs", "FreqHz", "ZData", "ZFit"
          }) &&
       Length[spec["Tau"]] == Length[spec["g"]] &&
       Length[spec["Tau"]] > 0 &&
       Length[spec["FreqHz"]] == Length[spec["ZData"]] &&
       Length[spec["FreqHz"]] == Length[spec["ZFit"]] &&
       VectorQ[N[spec["Tau"]], NumericQ] &&
       VectorQ[N[spec["g"]], NumericQ] &&
       NumericQ[N[spec["C0"]]] &&
       NumericQ[N[spec["Rs"]]];

     If[! ok,
      msg =
       If[
        spec === $Failed,
        "DCTSpectrum returned $Failed",
        If[
         AssociationQ[spec],
         "DCTSpectrum returned Association but failed validation. Keys = " <>
          ToString[Keys[spec]],
         "DCTSpectrum did not return an Association. Head = " <>
          ToString[Head[spec]]
        ]
       ];
      Print["Fit rejected for ", FileNameTake[file], ": ", msg];
     ];

     <|"File" -> file, "Spec" -> spec, "OK" -> ok, "Message" -> msg|>
    ],
    {f, txtFiles}
   ],
   Column[{
     Style["DCT batch progress", 14, Bold],
     Row[{
       ProgressIndicator[i/Max[nFiles, 1], {0, 1}],
       "  ",
       NumberForm[100. i/Max[nFiles, 1], {3, 1}],
       "%   (", i, "/", nFiles, ")"
      }],
     Row[{"Current file:  ", Style[currentFile, 12]}]
    }]
  ];

nOK = Count[results[[All, "OK"]], True];
Print["Succeeded on ", nOK, " / ", Length[txtFiles], " files."];

(* ============================================================ *)
(* 3) PLOT DCT RESULTS                                          *)
(* ============================================================ *)

goodSpecs = Select[results, #OK === True &];

If[Length[goodSpecs] == 0,
  Print["No successful fits to plot."];
  failures = Select[results, #OK =!= True &];
  If[Length[failures] > 0,
   Print["\nFailures:"];
   Do[
    Print["  ", FileNameTake[fail["File"]], " : ", fail["Message"]],
    {fail, failures}
   ];
  ];
  Abort[];
];

traces = (Transpose[{#Spec["Tau"], #Spec["g"]}] &) /@ goodSpecs;
labels = FileNameTake /@ (goodSpecs[[All, "File"]]);

colors =
  If[
   Length[traces] <= 1,
   {ColorData["Rainbow"][0.2]},
   Table[
    ColorData["Rainbow"][u],
    {u, 0.05, 0.95, (0.95 - 0.05)/(Length[traces] - 1)}
   ]
  ];

ListLinePlot[
 traces,
 Joined -> True,
 PlotStyle -> (Directive[#, AbsoluteThickness[2.2]] & /@ colors),
 ScalingFunctions -> {"Log10", None},
 Frame -> True,
 Axes -> False,
 FrameLabel -> {"\[Tau] (s)", "g(\[Tau]) (F/decade)"},
 PlotRange -> {Automatic, All},
 PlotLegends -> Placed[LineLegend[colors, labels], Right],
 ImageSize -> 700
]

(* ============================================================ *)
(* OPTIONAL: PRINT FAILURES                                     *)
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
(* Convert spectra to k-space                                   *)
(* ============================================================ *)

tracesK = (Transpose[{1/#Spec["Tau"], 10^6*#Spec["g"]}] &) /@ goodSpecs;

(* ============================================================ *)
(* Integrate g(k) over log10(k) to obtain distributed capacitance *)
(* ============================================================ *)

capFromGK[data_] := Module[
  {sorted, k, g, x},
  
  sorted = SortBy[data, First];
  k = sorted[[All, 1]];
  g = sorted[[All, 2]];
  
  x = Log10[k];
  
  If[Length[k] < 2, Return[0.0]];
  
  Total[
    Table[
      0.5*(g[[j]] + g[[j + 1]])*(x[[j + 1]] - x[[j]]),
      {j, 1, Length[k] - 1}
    ]
  ]
];

(* note: divide by 10^6 again since tracesK was scaled *)
cDistList = capFromGK /@ (Transpose[{#[[All,1]], #[[All,2]]/10^6}]& /@ tracesK);

(* ============================================================ *)
(* Extract C0 from fits                                         *)
(* ============================================================ *)

c0List = (#Spec["C0"] &) /@ goodSpecs;

cTotalList = cDistList + c0List;

(* ============================================================ *)
(* Plot styling                                                 *)
(* ============================================================ *)

n = Length[tracesK];

plotColors =
 Table[
  ColorData["Rainbow"][u],
  {u, 0.05, 0.95, (0.95 - 0.05)/Max[n - 1, 1]}
 ];

lineStyles = Map[Directive[#, AbsoluteThickness[2.5]] &, plotColors];

legendLabels =
 If[ValueQ[labels] && Length[labels] == n,
  labels,
  Table["Trace " <> ToString[i], {i, n}]
 ];

(* ============================================================ *)
(* Publication-quality DCT plot                                 *)
(* ============================================================ *)

dctPlot =
 ListLinePlot[
  tracesK,
  
  PlotStyle -> lineStyles,
  
  ScalingFunctions -> {"Log10", None},
  
  Frame -> True,
  Axes -> False,
  
  Background -> White,
  
  FrameStyle -> Directive[Black, Thickness[0.002]],
  
  FrameLabel -> {
    Style["Electron-transfer rate k (s^-1)", 17],
    Style["g(k) (\[Micro]F/decade)", 17]
  },
  
  LabelStyle -> Directive[Black, 16],
  
  PlotRange -> {Automatic, {0, 300}},
  
  ImageSize -> 800,
  
  PlotLegends -> Placed[
    LineLegend[lineStyles, legendLabels],
    Right
  ]
 ];

dctPlot

(* ============================================================ *)
(* Capacitance summary table                                    *)
(* ============================================================ *)

capSummary =
 Table[
  {
   legendLabels[[i]],
   10^6 * cDistList[[i]],
   10^6 * c0List[[i]],
   10^6 * cTotalList[[i]]
   },
  {i, 1, n}
 ];



(* ============================================================ *)
(* Filtered DCT plot: E0 only, colored by temperature            *)
(* ============================================================ *)

tracesKAll =
  (Transpose[{1/#Spec["Tau"], 10^6*#Spec["g"]}] &) /@ goodSpecs;

labelsAll =
  FileNameTake /@ (goodSpecs[[All, "File"]]);

e0Mask =
  StringContainsQ[#, "E0_"] & /@ labelsAll;

tracesKE0 = Pick[tracesKAll, e0Mask];
labelsE0 = Pick[labelsAll, e0Mask]

If[Length[tracesKE0] == 0,
  Print["No E0 traces found."];
  Abort[];
];

tempColor[label_] :=
  Which[
    StringContainsQ[label, "20C"], Yellow,
    StringContainsQ[label, "30C"], Orange,
    StringContainsQ[label, "40C"], Red,
    True, Gray
  ];
lineStylesE0 =
  Table[
    Which[
      StringContainsQ[labelsE0[[i]], "_1_"],
        Directive[tempColor[labelsE0[[i]]], AbsoluteThickness[2.5]],

      StringContainsQ[labelsE0[[i]], "_2_"],
        Directive[tempColor[labelsE0[[i]]], AbsoluteThickness[2.5], Dashed],

      StringContainsQ[labelsE0[[i]], "_3_"],
        Directive[tempColor[labelsE0[[i]]], AbsoluteThickness[2.5], DotDashed],

      True,
        Directive[tempColor[labelsE0[[i]]], AbsoluteThickness[2.5]]
    ],
    {i, Length[labelsE0]}
  ];
  
dctPlotE0 =
  ListLinePlot[
    tracesKE0,
    PlotStyle -> lineStylesE0,
    ScalingFunctions -> {"Log10", None},
    Frame -> True,
    Axes -> False,
    Background -> White,
    FrameStyle -> Directive[Black, Thickness[0.002]],
    FrameLabel -> {
      Style["Electron-transfer rate k (s^-1)", 17],
      Style["g(k) (\[Micro]F/decade)", 17]
    },
    LabelStyle -> Directive[Black, 16],
    PlotRange -> {Automatic, {0, 10}},
    ImageSize -> 800,
    PlotLegends -> Placed[
      LineLegend[lineStylesE0, labelsE0],
      Right
    ]
  ];

dctPlotE0


(* ============================================================ *)
(* E0 DCT plots: one panel per electrode, colored by temperature *)
(* ============================================================ *)

tracesKAll =
  (Transpose[{1/#Spec["Tau"], 10^6*#Spec["g"]}] &) /@ goodSpecs;

labelsAll = FileNameTake /@ (goodSpecs[[All, "File"]]);

e0Mask = StringContainsQ[#, "E0_"] & /@ labelsAll;

tracesKE0 = Pick[tracesKAll, e0Mask];
labelsE0 = Pick[labelsAll, e0Mask];

If[Length[tracesKE0] == 0,
  Print["No E0 traces found."];
  Abort[];
];

tempColor[label_] :=
  Which[
    StringContainsQ[label, "20C"], Yellow,
    StringContainsQ[label, "30C"], Orange,
    StringContainsQ[label, "40C"], Red,
    True, Gray
  ];

makeElectrodePlot[eTag_] := Module[
  {mask, traces, labels, styles},
  
  mask = StringContainsQ[#, eTag] & /@ labelsE0;
  traces = Pick[tracesKE0, mask];
  labels = Pick[labelsE0, mask];
  styles = Directive[tempColor[#], AbsoluteThickness[3]] & /@ labels;
  
  ListLinePlot[
    traces,
    PlotStyle -> styles,
    ScalingFunctions -> {"Log10", None},
    Frame -> True,
    Axes -> False,
    Background -> White,
    FrameStyle -> Directive[Black, Thickness[0.002]],
    FrameLabel -> {
      Style["k (s^-1)", 15],
      Style["g(k) (\[Mu]F/decade)", 15]
    },
    LabelStyle -> Directive[Black, 14],
    PlotRange -> {Automatic, {0, 10}},
    ImageSize -> 420,
    PlotLabel -> Style["Electrode " <> StringReplace[eTag, {"_" -> ""}], 16],
    PlotLegends -> Placed[
      LineLegend[styles, labels],
      Below
    ]
  ]
];

GraphicsGrid[
  {{
    makeElectrodePlot["_1_"],
    makeElectrodePlot["_2_"],
    makeElectrodePlot["_3_"]
  }},
  ImageSize -> 1300,
  Spacings -> {0.5, 0.5}
]


(* ============================================================ *)
(* Extract solution resistance from fits                        *)
(* ============================================================ *)

rsList = (#Spec["Rs"] &) /@ goodSpecs;

(* ============================================================ *)
(* Capacitance summary table                                    *)
(* ============================================================ *)

capSummary =
 Table[
  {
   legendLabels[[i]],
   rsList[[i]],
   10^6 * cDistList[[i]],
   10^6 * c0List[[i]],
   10^6 * cTotalList[[i]]
   },
  {i, 1, n}
 ];

Grid[
 Prepend[
  capSummary,
  {
   "Trace",
   "Rs (\[CapitalOmega])",
   "\[Integral] g(k) dlog10k (\[Micro]F)",
   "C0 (\[Micro]F)",
   "Total C (\[Micro]F)"
   }
  ],
 Frame -> All,
 Alignment -> Left,
 ItemStyle -> Directive[Black, 13]
]


(* ============================================================ *)
(* Convert spectra to k-space                                   *)
(* ============================================================ *)

tracesK = (Transpose[{1/#Spec["Tau"], #Spec["g"]}] &) /@ goodSpecs;

(* ============================================================ *)
(* Build cumulative capacitance vs k from g(k)                  *)
(* cumulative C(k) = Integral g(k') dlog10(k') from low k up    *)
(* ============================================================ *)

cumCapVsK[data_] := Module[
  {sorted, k, g, x, c},
  
  sorted = SortBy[data, First];
  k = sorted[[All, 1]];
  g = sorted[[All, 2]];
  
  If[Length[k] < 2 || Length[g] != Length[k], Return[{}]];
  
  x = Log10[k];
  
  c = Prepend[
    Accumulate[
      Table[
        0.5*(g[[j]] + g[[j + 1]])*(x[[j + 1]] - x[[j]]),
        {j, 1, Length[k] - 1}
      ]
    ],
    0
  ];
  
  Transpose[{k, 10^6*c}]
];

cumTracesK = cumCapVsK /@ tracesK;

(* ============================================================ *)
(* Extract C0 from fits                                         *)
(* ============================================================ *)

c0List = (#Spec["C0"] &) /@ goodSpecs;

(* ============================================================ *)
(* Plot styling                                                 *)
(* ============================================================ *)

n = Length[cumTracesK];

plotColors =
 Table[
  ColorData["Rainbow"][u],
  {u, 0.05, 0.95, (0.95 - 0.05)/Max[n - 1, 1]}
 ];

lineStyles = Map[Directive[#, AbsoluteThickness[2.5]] &, plotColors];

legendLabels =
 If[ValueQ[labels] && Length[labels] == n,
  labels,
  Table["Trace " <> ToString[i], {i, n}]
 ];

(* ============================================================ *)
(* Publication-quality cumulative capacitance plot              *)
(* ============================================================ *)

cumPlot =
 ListLinePlot[
  cumTracesK,
  
  PlotStyle -> lineStyles,
  
  ScalingFunctions -> {"Log10", None},
  
  Frame -> True,
  Axes -> False,
  
  Background -> White,
  
  FrameStyle -> Directive[Black, Thickness[0.002]],
  
  FrameLabel -> {
    Style["Electron-transfer rate k (s^-1)", 17],
    Style["Cumulative capacitance (\[Micro]F)", 17]
  },
  
  LabelStyle -> Directive[Black, 16],
  
  PlotRange -> All,
  
  ImageSize -> 800,
  
  PlotLegends -> Placed[
    LineLegend[lineStyles, legendLabels],
    Right
  ]
 ];

cumPlot

(* ============================================================ *)
(* Capacitance summary table                                    *)
(* ============================================================ *)

cDistList = Module[{x, y},
    x = #[[All, 1]];
    y = #[[All, 2]];
    If[Length[y] > 0, Last[y]/10^6, 0.]
  ] & /@ cumTracesK;

cTotalList = cDistList + c0List;

capSummary =
 Table[
  {
   legendLabels[[i]],
   10^6*cDistList[[i]],
   10^6*c0List[[i]],
   10^6*cTotalList[[i]]
   },
  {i, 1, n}
 ];

Grid[
 Prepend[
  capSummary,
  {
   "Trace",
   "\[Integral] g(k) dlog10k (\[Micro]F)",
   "C0 (\[Micro]F)",
   "Total C (\[Micro]F)"
   }
  ],
 Frame -> All,
 Alignment -> Left,
 ItemStyle -> Directive[Black, 13]
]


(#Spec["CPEAlpha"] &) /@ goodSpecs
(#Spec["Y0"] &) /@ goodSpecs
(#Spec["C0"] &) /@ goodSpecs
(#Spec["Rs"] &) /@ goodSpecs
Min[(1/#Spec["Tau"] &) /@ goodSpecs]
Max[(1/#Spec["Tau"] &) /@ goodSpecs]
Min[(1/#Spec["Tau"] &) /@ goodSpecs]/(2 \[Pi])
Max[(1/#Spec["Tau"] &) /@ goodSpecs]/(2 \[Pi])


(* helper function *)
getElectrodeFromFile[file_String] := Module[{name, hit},
  name = FileNameTake[file];
  hit = StringCases[
    name,
    RegularExpression["E(\\d+)"] :> ("E" <> "$1")
  ];
  If[hit === {}, Missing["NoElectrode"], First[hit]]
];


goodSpecsE1 =
  Select[
    goodSpecs,
    getElectrodeFromFile[#["File"]] === "E1" &
  ];

goodSpecsE1 = goodSpecs;

(* elapsed time in hours *)
expTimesHr = (#Spec["FinishTimeS"]/3600) & /@ goodSpecsE1;
elapsedHr = expTimesHr - Min[expTimesHr];

(* choose target legend times *)
targetTimesHr = Range[0, 24, 4];

(* for each target time, pick the closest dataset *)
closestIdx =
  DeleteDuplicates[
    (First @ Ordering[Abs[elapsedHr - #], 1]) & /@ targetTimesHr
  ];

goodSpecsE1 = goodSpecsE1[[closestIdx]];
elapsedHrSelected = elapsedHr[[closestIdx]];

labelsE1 = (ToString[#] <> " hrs") & /@ targetTimesHr[[;; Length[closestIdx]]];

tracesK =
  (Transpose[{1/#Spec["Tau"], #Spec["g"]}] &) /@
    goodSpecsE1;

nE1 = Length[tracesK];

rainbowColors =
  If[
    nE1 <= 1,
    {ColorData["Rainbow"][0.2]},
    Table[
      ColorData["Rainbow"][u],
      {u, 0.1, 0.9, (0.9 - 0.1)/(nE1 - 1)}
    ]
  ];

legend =
  Framed[
    LineLegend[
      Map[Directive[#, AbsoluteThickness[2.2]] &, rainbowColors],
      labelsE1,
      LegendMarkerSize -> 28,
      LabelStyle -> Directive[Black, 11, FontFamily -> "Arial"],
      LegendLayout -> "Column",
      Spacings -> 0.2
    ],
    Background -> White,
    FrameStyle -> White,
    FrameMargins -> 6
  ];

(* ----- powers-of-10 x ticks ----- *)
kMin = Min[Flatten[(1/#Spec["Tau"]) & /@ goodSpecsE1]];
kMax = Max[Flatten[(1/#Spec["Tau"]) & /@ goodSpecsE1]];

xTicks =
  Table[
    {10.^p, Superscript[10, p]},
    {p, Floor[Log10[kMin]], Ceiling[Log10[kMax]]}
  ];

plotE1 =
  Legended[
    ListLinePlot[
      tracesK,
      PlotStyle -> Map[Directive[#, AbsoluteThickness[2.2]] &, rainbowColors],
      ScalingFunctions -> {"Log10", None},
      Frame -> True,
      Axes -> False,
      FrameTicks -> {{Automatic, None}, {xTicks, None}},
      Background -> White,
      FrameStyle -> Directive[Black, AbsoluteThickness[1.2]],
      FrameTicksStyle -> Directive[Black, 18],
      LabelStyle -> Directive[Black, 18, FontFamily -> "Arial"],
      FrameLabel -> {
        Style[
          Row[{"Electron transfer rate, k (", Superscript["s", -1], ")"}],
          18, Black, FontFamily -> "Arial"
        ],
        Style["g(k) (F/decade)", 18, Black, FontFamily -> "Arial"]
      },
      PlotRange -> {Automatic, {0, 300 10^-6}},
      PlotRangePadding -> {{Scaled[0.02], Scaled[0.02]}, {Scaled[0.02], Scaled[0.04]}},
      ImageSize -> 400,
      AspectRatio -> 0.68,
      BaseStyle -> {FontFamily -> "Arial", 16},
      GridLines -> None
    ],
    Placed[legend, Right]
  ];

plotE1


Export["driftingPeaks.png",plotE1]


  



goodSpecsE1 = Select[ goodSpecs, getElectrodeFromFile[#["File"]] === "E1" & ];
rawLabelsE1 = FileNameTake /@ (goodSpecsE1[[All, "File"]]);

keepMask = Map[
  StringContainsQ[#, "_0uM"] ||
  StringContainsQ[#, "_100uM"] ||
  StringContainsQ[#, "_250uM"] ||
  StringContainsQ[#, "_562uM"] &,
  rawLabelsE1
];

goodSpecsE1 = Pick[goodSpecsE1, keepMask];
rawLabelsE1 = Pick[rawLabelsE1, keepMask];

(* extract concentration values *)
concVals =
  ToExpression @ StringReplace[
    StringCases[#, DigitCharacter .. ~~ "uM"][[1]],
    "uM" -> ""
  ] & /@ rawLabelsE1;

(* sort by concentration descending *)
sortIdx = Reverse @ Ordering[concVals];

goodSpecsE1 = goodSpecsE1[[sortIdx]];
rawLabelsE1 = rawLabelsE1[[sortIdx]];
concVals = concVals[[sortIdx]];

labelsE1 = (ToString[#] <> " \[Micro]M") & /@ concVals;

tracesK =
  (Transpose[{1/#Spec["Tau"], #Spec["g"]}] &) /@
    goodSpecsE1;

nE1 = Length[tracesK];

rainbowColors =
  If[
    nE1 <= 1,
    {ColorData["DeepSeaColors"][0.2]},
    Table[
      ColorData["DeepSeaColors"][u],
      {u, 0.1, 0.9, (0.9 - 0.1)/(nE1 - 1)}
    ]
  ];

legend =
  Framed[
    LineLegend[
      Map[Directive[#, AbsoluteThickness[2.2]] &, rainbowColors],
      labelsE1,
      LegendMarkerSize -> 28,
      LabelStyle -> Directive[Black, 11, FontFamily -> "Arial"],
      LegendLayout -> "Column"
    ],
    Background -> White,
    FrameStyle -> White,
    FrameMargins -> 6
  ];

(* ----- powers-of-10 x ticks ----- *)
kMin = Min[Flatten[(1/#Spec["Tau"]) & /@ goodSpecsE1]];
kMax = Max[Flatten[(1/#Spec["Tau"]) & /@ goodSpecsE1]];

xTicks =
  Table[
    {10.^p, Superscript[10, p]},
    {p, Floor[Log10[kMin]], Ceiling[Log10[kMax]]}
  ];

plotE1 =
  Legended[
    ListLinePlot[
      tracesK,
      PlotStyle -> Map[Directive[#, AbsoluteThickness[2.2]] &, rainbowColors],
      ScalingFunctions -> {"Log10", None},
      Frame -> True,
      Axes -> False,
      FrameTicks -> {{Automatic, None}, {xTicks, None}},
      Background -> White,
      FrameStyle -> Directive[Black, AbsoluteThickness[1.2]],
      FrameTicksStyle -> Directive[Black, 18],
      LabelStyle -> Directive[Black, 18, FontFamily -> "Arial"],
      FrameLabel -> {
        Style[
          Row[{"Electron transfer rate, k (", Superscript["s", -1], ")"}],
          18, Black, FontFamily -> "Arial"
        ],
        Style["g(k) (F/decade)", 18, Black, FontFamily -> "Arial"]
      },
      PlotRange -> {Automatic, {0,1.*10^-6}},
      PlotRangePadding -> {{Scaled[0.02], Scaled[0.02]}, {Scaled[0.02], Scaled[0.04]}},
      ImageSize -> 600,
      AspectRatio -> 0.68,
      BaseStyle -> {FontFamily -> "Arial", 16},
      GridLines -> None
    ],
    Placed[legend, {0.85, 0.5}]
  ];

plotE1


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


Export["titrationPeaks.png",plotE1]


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
(* Rainbow Bode overlays: RAW data for ALL files                *)
(*   Col 1: RAW                                                 *)
(*   Col 2: FIT                                                 *)
(* ============================================================ *)
(* 
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
*)



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
(* ---------- RAW impedance ---------- *)
rawZResults =
	Table[
		Module[
			{
				file, raw, fRaw, zRaw,
				magTrace, phaseTrace, reTrace, imTrace, ok, msg
			},

			file = files[[k]];
			raw  = Quiet @ Check[DCT`Private`importEISTxt[file], $Failed];

			If[
				raw === $Failed || FailureQ[raw] || !AssociationQ[raw] ||
				!KeyExistsQ[raw, "FreqHz"] || !KeyExistsQ[raw, "Z"],
				<|
					"File" -> file,
					"Label" -> labels[[k]],
					"OK" -> False,
					"Message" -> "raw parse failed"
				|>,
				
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
			]
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
			raw  = Quiet @ Check[DCT`Private`importEISTxt[file], $Failed];

			If[
				raw === $Failed || FailureQ[raw] || !AssociationQ[raw] || !KeyExistsQ[raw, "FreqHz"],
				<|
					"File" -> file,
					"Label" -> labels[[k]],
					"OK" -> False,
					"Message" -> "raw parse failed or missing FreqHz"
				|>,

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
					<|
						"File" -> file,
						"Label" -> labels[[k]],
						"OK" -> False,
						"Message" -> "missing fit parameters for impedance reconstruction"
					|>,

					If[Min[Abs[YintFit]] < 10^-30,
						<|
							"File" -> file,
							"Label" -> labels[[k]],
							"OK" -> False,
							"Message" -> "YintFit too small"
						|>,

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
					]
				]
			]
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
(* Admittance Bode overlays for TOTAL SYSTEM                    *)
(*   Col 1: RAW total admittance                                *)
(*   Col 2: FIT total admittance                                *)
(*   Col 3: ABS residual error                                  *)
(*   Row 1: |Y|                                                 *)
(*   Row 2: Phase(Y)                                            *)
(*   Row 3: Re(Y)                                               *)
(*   Row 4: -Im(Y)                                              *)
(* accommodates useConstantPhaseElementFlag                     *)
(* ============================================================ *)

If[Length[goodSpecs] == 0,
	Print["No goodSpecs to plot."];
	Abort[];
];

labels = FileNameTake /@ (goodSpecs[[All, "File"]]);
files  = goodSpecs[[All, "File"]];

(* ---------- RAW total admittance ---------- *)
rawYResults =
	Table[
		Module[
			{
				file, raw, fRaw, zRaw, yRaw,
				magTrace, phaseTrace, reTrace, imTrace, ok, msg
			},

			file = files[[k]];
			raw  = Quiet @ Check[DCT`Private`importEISTxt[file], $Failed];

			If[
				raw === $Failed || FailureQ[raw] || !AssociationQ[raw] ||
				!KeyExistsQ[raw, "FreqHz"] || !KeyExistsQ[raw, "Z"],
				<|
					"File" -> file,
					"Label" -> labels[[k]],
					"OK" -> False,
					"Message" -> "raw parse failed"
				|>,

				fRaw = raw["FreqHz"];
				zRaw = raw["Z"];

				If[Min[Abs[zRaw]] < 10^-15,
					<|
						"File" -> file,
						"Label" -> labels[[k]],
						"OK" -> False,
						"Message" -> "Zraw too small"
					|>,

					yRaw = 1/zRaw;

					magTrace   = Transpose[{fRaw, Abs[yRaw]}];
					phaseTrace = Transpose[{fRaw, (180./Pi) * Arg[yRaw]}];
					reTrace    = Transpose[{fRaw, Re[yRaw]}];
					imTrace    = Transpose[{fRaw, -Im[yRaw]}];

					ok = VectorQ[magTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &] &&
						 VectorQ[phaseTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &] &&
						 VectorQ[reTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &] &&
						 VectorQ[imTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &];

					msg = If[ok, "", "raw admittance trace non-numeric"];

					<|
						"File" -> file,
						"Label" -> labels[[k]],
						"FreqHz" -> fRaw,
						"YRaw" -> yRaw,
						"OK" -> ok,
						"MagTrace" -> magTrace,
						"PhaseTrace" -> phaseTrace,
						"ReTrace" -> reTrace,
						"ImTrace" -> imTrace,
						"Message" -> msg
					|>
				]
			]
		],
		{k, Length[files]}
	];

goodRawY = Select[rawYResults, TrueQ[#["OK"]] &];
badRawY  = Select[rawYResults, !TrueQ[#["OK"]] &];

Print["Raw total admittance overlay: parsed OK = ", Length[goodRawY], " / ", Length[rawYResults]];

If[Length[badRawY] > 0,
	Print["Failures:"];
	Scan[
		(Print["\t", FileNameTake[#["File"]], " : ", #["Message"]]) &,
		badRawY
	];
];

If[Length[goodRawY] == 0,
	Print["No raw total admittance spectra could be parsed."];
	Abort[];
];

rawYMagTraces   = goodRawY[[All, "MagTrace"]];
rawYPhaseTraces = goodRawY[[All, "PhaseTrace"]];
rawYReTraces    = goodRawY[[All, "ReTrace"]];
rawYImTraces    = goodRawY[[All, "ImTrace"]];
rawYLabels      = goodRawY[[All, "Label"]];
rawYColors      = ColorData["Rainbow"] /@ Rescale[Range[Length[goodRawY]]];

(* ---------- FIT total admittance ---------- *)
fullFitYResults =
	Table[
		Module[
			{
				spec, raw, tau, g, c0, y0, alphaCPE, useCPE, rs,
				fRaw, \[Omega], \[CapitalDelta]log, Kmat, YintFit, Zfit, Yfit,
				magTrace, phaseTrace, reTrace, imTrace, file, ok, msg
			},

			spec = goodSpecs[[k, "Spec"]];
			file = goodSpecs[[k, "File"]];
			raw  = Quiet @ Check[DCT`Private`importEISTxt[file], $Failed];

			If[
				raw === $Failed || FailureQ[raw] || !AssociationQ[raw] || !KeyExistsQ[raw, "FreqHz"],
				<|
					"File" -> file,
					"Label" -> labels[[k]],
					"OK" -> False,
					"Message" -> "raw parse failed or missing FreqHz"
				|>,

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
					<|
						"File" -> file,
						"Label" -> labels[[k]],
						"OK" -> False,
						"Message" -> "missing fit parameters for admittance reconstruction"
					|>,

					If[Min[Abs[YintFit]] < 10^-30,
						<|
							"File" -> file,
							"Label" -> labels[[k]],
							"OK" -> False,
							"Message" -> "YintFit too small"
						|>,

						Zfit = rs + 1/YintFit;

						If[Min[Abs[Zfit]] < 10^-15,
							<|
								"File" -> file,
								"Label" -> labels[[k]],
								"OK" -> False,
								"Message" -> "Zfit too small"
							|>,

							Yfit = 1/Zfit;

							magTrace   = Transpose[{fRaw, Abs[Yfit]}];
							phaseTrace = Transpose[{fRaw, (180./Pi) * Arg[Yfit]}];
							reTrace    = Transpose[{fRaw, Re[Yfit]}];
							imTrace    = Transpose[{fRaw, -Im[Yfit]}];

							ok = VectorQ[magTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &] &&
								 VectorQ[phaseTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &] &&
								 VectorQ[reTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &] &&
								 VectorQ[imTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &];

							msg = If[ok, "", "fit admittance trace non-numeric"];

							<|
								"File" -> file,
								"Label" -> labels[[k]],
								"FreqHz" -> fRaw,
								"YFit" -> Yfit,
								"MagTrace" -> magTrace,
								"PhaseTrace" -> phaseTrace,
								"ReTrace" -> reTrace,
								"ImTrace" -> imTrace,
								"OK" -> ok,
								"Message" -> msg
							|>
						]
					]
				]
			]
		],
		{k, Length[goodSpecs]}
	];

goodFullFitY = Select[fullFitYResults, TrueQ[#["OK"]] &];
badFullFitY  = Select[fullFitYResults, !TrueQ[#["OK"]] &];

If[Length[badFullFitY] > 0,
	Print["Failures:"];
	Scan[
		(Print["\t", #["File"], " : ", #["Message"]]) &,
		badFullFitY
	];
];

If[Length[goodFullFitY] == 0,
	Print["No fit total admittance spectra could be computed."];
	Abort[];
];

fullFitYMagTraces   = goodFullFitY[[All, "MagTrace"]];
fullFitYPhaseTraces = goodFullFitY[[All, "PhaseTrace"]];
fullFitYReTraces    = goodFullFitY[[All, "ReTrace"]];
fullFitYImTraces    = goodFullFitY[[All, "ImTrace"]];
fitYColors          = ColorData["Rainbow"] /@ Rescale[Range[Length[goodFullFitY]]];

(* ---------- Residual traces ---------- *)
residualYResults =
	Table[
		Module[
			{
				rawAssoc, fitAssoc, fRaw, yRaw, yFit,
				magResTrace, phaseResTrace, reResTrace, imResTrace
			},

			rawAssoc = rawYResults[[k]] /. HoldPattern[Return[x_]] :> x;
			fitAssoc = fullFitYResults[[k]] /. HoldPattern[Return[x_]] :> x;

			If[!TrueQ[Lookup[rawAssoc, "OK", False]] || !TrueQ[Lookup[fitAssoc, "OK", False]],
				<|
					"OK" -> False,
					"Message" -> "raw/fit pair unavailable"
				|>,

				fRaw = rawAssoc["FreqHz"];
				yRaw = rawAssoc["YRaw"];
				yFit = fitAssoc["YFit"];

				If[Length[fRaw] =!= Length[yRaw] || Length[yRaw] =!= Length[yFit],
					<|
						"OK" -> False,
						"Message" -> "raw and fit lengths mismatch"
					|>,

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
				]
			]
		],
		{k, Length[files]}
	];

goodResidualY = Select[residualYResults, TrueQ[#["OK"]] &];

If[Length[goodResidualY] == 0,
	Print["No residual total admittance traces could be computed."];
	Abort[];
];

residualYMagTraces   = goodResidualY[[All, "MagResTrace"]];
residualYPhaseTraces = goodResidualY[[All, "PhaseResTrace"]];
residualYReTraces    = goodResidualY[[All, "ReResTrace"]];
residualYImTraces    = goodResidualY[[All, "ImResTrace"]];
resYColors           = ColorData["Rainbow"] /@ Rescale[Range[Length[goodResidualY]]];

(* ---------- Shared y-axis ranges ---------- *)

allMagY = Select[
	Join[
		Flatten[rawYMagTraces[[All, All, 2]]],
		Flatten[fullFitYMagTraces[[All, All, 2]]]
	],
	NumericQ[#] && # > 0 &
];

allPhaseY = Select[
	Join[
		Flatten[rawYPhaseTraces[[All, All, 2]]],
		Flatten[fullFitYPhaseTraces[[All, All, 2]]]
	],
	NumericQ
];

allReY = Select[
	Join[
		Flatten[rawYReTraces[[All, All, 2]]],
		Flatten[fullFitYReTraces[[All, All, 2]]]
	],
	NumericQ
];

allImY = Select[
	Join[
		Flatten[rawYImTraces[[All, All, 2]]],
		Flatten[fullFitYImTraces[[All, All, 2]]]
	],
	NumericQ
];

allMagResY = Select[Flatten[residualYMagTraces[[All, All, 2]]], NumericQ[#] && # > 0 &];
allPhaseResY = Select[Flatten[residualYPhaseTraces[[All, All, 2]]], NumericQ];
allReResY = Select[Flatten[residualYReTraces[[All, All, 2]]], NumericQ];
allImResY = Select[Flatten[residualYImTraces[[All, All, 2]]], NumericQ];

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
phaseResYRange = {Min[allPhaseResY], Min[10,Max[allPhaseResY]]};
reResYRange    = {Min[allReResY], Max[allReResY]};
imResYRange    = {Min[allImResY], Max[allImResY]};

Print["Shared |Y| y-range = ", magYRange];
Print["Shared phase y-range = ", phaseYRange];
Print["Shared Re(Y) y-range = ", reYRange];
Print["Shared -Im(Y) y-range = ", imYRange];

(* ---------- Plots ---------- *)

rawYMagPlot =
	ListLogLogPlot[
		rawYMagTraces,
		PlotStyle -> rawYColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "|Yraw| (S)"},
		PlotRange -> {All, magYRange},
		ImageSize -> 900,
		Background -> White,
		PlotLabel -> Row[{"RAW |Y|   n=", Length[goodRawY]}]
	];

fitYMagPlot =
	ListLogLogPlot[
		fullFitYMagTraces,
		PlotStyle -> fitYColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "|Yfit| (S)"},
		PlotRange -> {All, magYRange},
		ImageSize -> 900,
		Background -> White,
		PlotLabel -> Row[{"FIT |Y|   n=", Length[goodFullFitY]}]
	];

resYMagPlot =
	ListLogLogPlot[
		residualYMagTraces,
		PlotStyle -> resYColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "|(Yraw - Yfit)/Yraw|"},
		PlotRange -> {All, magResYRange},
		ImageSize -> 900,
		Background -> White,
		PlotLabel -> "ABS residual |(Yraw - Yfit)/Yraw|"
	];

rawYPhasePlot =
	ListLogLinearPlot[
		rawYPhaseTraces,
		PlotStyle -> rawYColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "Phase(Yraw) (deg)"},
		PlotRange -> {All, phaseYRange},
		ImageSize -> 900,
		Background -> White,
		PlotLabel -> Row[{"RAW phase(Y)   n=", Length[goodRawY]}]
	];

fitYPhasePlot =
	ListLogLinearPlot[
		fullFitYPhaseTraces,
		PlotStyle -> fitYColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "Phase(Yfit) (deg)"},
		PlotRange -> {All, phaseYRange},
		ImageSize -> 900,
		Background -> White,
		PlotLabel -> Row[{"FIT phase(Y)   n=", Length[goodFullFitY]}]
	];

resYPhasePlot =
	ListLogLinearPlot[
		residualYPhaseTraces,
		PlotStyle -> resYColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "|Phase raw - fit| (deg)"},
		PlotRange -> {All, phaseResYRange},
		ImageSize -> 900,
		Background -> White,
		PlotLabel -> "ABS phase residual"
	];

rawYRePlot =
	ListLogLinearPlot[
		rawYReTraces,
		PlotStyle -> rawYColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "Re(Yraw) (S)"},
		PlotRange -> {All, reYRange},
		ImageSize -> 900,
		Background -> White,
		PlotLabel -> "RAW Re(Y)"
	];

fitYRePlot =
	ListLogLinearPlot[
		fullFitYReTraces,
		PlotStyle -> fitYColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "Re(Yfit) (S)"},
		PlotRange -> {All, reYRange},
		ImageSize -> 900,
		Background -> White,
		PlotLabel -> "FIT Re(Y)"
	];

resYRePlot =
	ListLogLinearPlot[
		residualYReTraces,
		PlotStyle -> resYColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "|Re raw - fit| (S)"},
		PlotRange -> {All, reResYRange},
		ImageSize -> 900,
		Background -> White,
		PlotLabel -> "ABS Re residual"
	];

rawYImPlot =
	ListLogLinearPlot[
		rawYImTraces,
		PlotStyle -> rawYColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "-Im(Yraw) (S)"},
		PlotRange -> {All, imYRange},
		ImageSize -> 900,
		Background -> White,
		PlotLabel -> "RAW -Im(Y)"
	];

fitYImPlot =
	ListLogLinearPlot[
		fullFitYImTraces,
		PlotStyle -> fitYColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "-Im(Yfit) (S)"},
		PlotRange -> {All, imYRange},
		ImageSize -> 900,
		Background -> White,
		PlotLabel -> "FIT -Im(Y)"
	];

resYImPlot =
	ListLogLinearPlot[
		residualYImTraces,
		PlotStyle -> resYColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "|(-Im raw) - (-Im fit)| (S)"},
		PlotRange -> {All, imResYRange},
		ImageSize -> 900,
		Background -> White,
		PlotLabel -> "ABS -Im residual"
	];

GraphicsGrid[
	{
		{rawYMagPlot, fitYMagPlot, resYMagPlot},
		{rawYPhasePlot, fitYPhasePlot, resYPhasePlot}
		,{rawYRePlot, fitYRePlot, resYRePlot},
		{rawYImPlot, fitYImPlot, resYImPlot}
	},
	ImageSize -> 1200
]


(* ============================================================ *)
(* Ratio diagnostic in TOTAL admittance space                   *)
(*   Magnitude: |Yraw / Yfit|                                   *)
(*   Phase:     Phase(Yraw) - Phase(Yfit)                       *)
(* ============================================================ *)

If[Length[rawYResults] == 0 || Length[fullFitYResults] == 0,
	Print["Need rawYResults and fullFitYResults first."];
	Abort[];
];

rawYResultsClean     = rawYResults /. HoldPattern[Return[x_]] :> x;
fullFitYResultsClean = fullFitYResults /. HoldPattern[Return[x_]] :> x;

ratioYResults =
	Table[
		Module[
			{
				rawAssoc, fitAssoc, fRaw, yRaw, yFit, yRatio,
				magTrace, phaseDiffTrace
			},

			rawAssoc = rawYResultsClean[[k]];
			fitAssoc = fullFitYResultsClean[[k]];

			If[
				!AssociationQ[rawAssoc] || !AssociationQ[fitAssoc],
				<|
					"OK" -> False,
					"Message" -> "raw/fit entry is not an Association"
				|>,

				If[
					!TrueQ[Lookup[rawAssoc, "OK", False]] || !TrueQ[Lookup[fitAssoc, "OK", False]],
					<|
						"OK" -> False,
						"Message" -> "raw/fit pair unavailable"
					|>,

					fRaw = Lookup[rawAssoc, "FreqHz", Missing["NoFreq"]];
					yRaw = Lookup[rawAssoc, "YRaw", Missing["NoYRaw"]];
					yFit = Lookup[fitAssoc, "YFit", Missing["NoYFit"]];

					If[
						MemberQ[{fRaw, yRaw, yFit}, _Missing],
						<|
							"OK" -> False,
							"Message" -> "missing FreqHz, YRaw, or YFit"
						|>,

						If[
							Length[fRaw] =!= Length[yRaw] || Length[yRaw] =!= Length[yFit],
							<|
								"OK" -> False,
								"Message" -> "raw and fit lengths mismatch"
							|>,

							If[
								Min[Abs[yFit]] < 10^-30,
								<|
									"OK" -> False,
									"Message" -> "Yfit too small for ratio"
								|>,

								yRatio = yRaw/yFit;

								magTrace       = Transpose[{fRaw, Abs[yRatio]}];
								phaseDiffTrace = Transpose[{fRaw, (180./Pi) * (Arg[yRaw] - Arg[yFit])}];

								<|
									"OK" -> True,
									"MagTrace" -> magTrace,
									"PhaseDiffTrace" -> phaseDiffTrace
								|>
							]
						]
					]
				]
			]
		],
		{k, Length[rawYResultsClean]}
	];

goodRatioY = Select[ratioYResults, TrueQ[Lookup[#, "OK", False]] &];
badRatioY  = Select[ratioYResults, !TrueQ[Lookup[#, "OK", False]] &];

If[Length[badRatioY] > 0,
	Print["Failures in total-Y ratio traces:"];
	MapIndexed[
		Function[{assoc, idx},
			Print["\ttrace ", First[idx], " : ", Lookup[assoc, "Message", "no message"]]
		],
		badRatioY
	];
];

If[Length[goodRatioY] == 0,
	Print["No total-Y ratio traces could be computed."];
	Abort[];
];

ratioYMagTraces    = goodRatioY[[All, "MagTrace"]];
phaseDiffYTraces   = goodRatioY[[All, "PhaseDiffTrace"]];
ratioYColors       = ColorData["Rainbow"] /@ Rescale[Range[Length[goodRatioY]]];

allRatioMagY = Select[
	Flatten[ratioYMagTraces[[All, All, 2]]],
	NumericQ[#] && # > 0 &
];

allPhaseDiffY = Select[
	Flatten[phaseDiffYTraces[[All, All, 2]]],
	NumericQ
];

If[allRatioMagY === {} || allPhaseDiffY === {},
	Print["Could not determine total-Y ratio plot ranges."];
	Abort[];
];

ratioMagYRange  = {Min[allRatioMagY], Max[allRatioMagY]};
phaseDiffYRange = {Min[allPhaseDiffY], Max[allPhaseDiffY]};

magRef = LogLinearPlot[
	1,
	{x,
	 Min[Flatten[ratioYMagTraces[[All, All, 1]]]],
	 Max[Flatten[ratioYMagTraces[[All, All, 1]]]]
	},
	PlotStyle -> Directive[Black, Dashed]
];

phaseRef = LogLinearPlot[
	0,
	{x,
	 Min[Flatten[phaseDiffYTraces[[All, All, 1]]]],
	 Max[Flatten[phaseDiffYTraces[[All, All, 1]]]]
	},
	PlotStyle -> Directive[Black, Dashed]
];

ratioYMagPlot =
	Show[
		ListLogLinearPlot[
			ratioYMagTraces,
			PlotStyle -> ratioYColors,
			Joined -> True,
			Frame -> True,
			Axes -> False,
			FrameLabel -> {"Frequency (Hz)", "|Yraw / Yfit|"},
			(* PlotRange -> {All, ratioMagYRange}, *)
			PlotRange -> {All, {0.97,1.03}},
			ImageSize -> 800,
			Background -> White,
			PlotLabel -> "Magnitude ratio for total-system admittance"
		],
		magRef
	];

phaseDiffYPlot =
	Show[
		ListLogLinearPlot[
			phaseDiffYTraces,
			PlotStyle -> ratioYColors,
			Joined -> True,
			Frame -> True,
			Axes -> False,
			FrameLabel -> {"Frequency (Hz)", "Phase(Yraw) - Phase(Yfit) (deg)"},
			(* PlotRange -> {All, phaseDiffYRange}, *)
			PlotRange -> {All, {-2,2.}},
			ImageSize -> 800,
			Background -> White,
			PlotLabel -> "Phase difference for total-system admittance"
		],
		phaseRef
	];

GraphicsRow[
	{ratioYMagPlot, phaseDiffYPlot},
	ImageSize -> 1200
]


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ============================================================ *)
(* RAW complex capacitance overlays                             *)
(* Ccomplex = 1/(I w Z)                                          *)
(*   Row 1: |Ccomplex|                                           *)
(*   Row 2: Phase(Ccomplex)                                      *)
(*   Row 3: Re(Ccomplex)                                         *)
(*   Row 4: Im(Ccomplex)                                         *)
(* ============================================================ *)

If[Length[goodSpecs] == 0,
	Print["No goodSpecs to plot."];
	Abort[];
];

labels = FileNameTake /@ (goodSpecs[[All, "File"]]);
files  = goodSpecs[[All, "File"]];
rawCResults =
	Table[
		Module[
			{
				file, raw, fRaw, zRaw, omega, cRaw,
				magTrace, phaseTrace, reTrace, imTrace, ok, msg
			},

			file = files[[k]];
			raw  = Quiet @ Check[DCT`Private`importEISTxt[file], $Failed];

			If[
				raw === $Failed || FailureQ[raw] || !AssociationQ[raw] ||
				!KeyExistsQ[raw, "FreqHz"] || !KeyExistsQ[raw, "Z"],
				<|
					"File" -> file,
					"Label" -> labels[[k]],
					"OK" -> False,
					"Message" -> "raw parse failed"
				|>,

				fRaw = raw["FreqHz"];
				zRaw = raw["Z"];
				omega = 2 Pi fRaw;

				If[
					Length[fRaw] =!= Length[zRaw] ||
					Min[Abs[omega]] <= 0 ||
					Min[Abs[zRaw]] <= 10^-30,
					<|
						"File" -> file,
						"Label" -> labels[[k]],
						"OK" -> False,
						"Message" -> "invalid raw data for capacitance conversion"
					|>,

					cRaw = 1/(I*omega*zRaw);

					magTrace   = Transpose[{fRaw, 10^6*Abs[cRaw]}];
					phaseTrace = Transpose[{fRaw, (180./Pi)*Arg[cRaw]}];
					reTrace    = Transpose[{fRaw, 10^6*Re[cRaw]}];
					imTrace    = Transpose[{fRaw, -10^6*Im[cRaw]}];

					ok = VectorQ[magTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &] &&
						 VectorQ[phaseTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &] &&
						 VectorQ[reTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &] &&
						 VectorQ[imTrace, MatchQ[#, {_?NumericQ, _?NumericQ}] &];

					msg = If[ok, "", "raw complex capacitance trace non-numeric"];

					<|
						"File" -> file,
						"Label" -> labels[[k]],
						"FreqHz" -> fRaw,
						"CRaw" -> cRaw,
						"OK" -> ok,
						"MagTrace" -> magTrace,
						"PhaseTrace" -> phaseTrace,
						"ReTrace" -> reTrace,
						"ImTrace" -> imTrace,
						"Message" -> msg
					|>
				]
			]
		],
		{k, Length[files]}
	];

goodRawC = Select[rawCResults, TrueQ[#["OK"]] &];
badRawC  = Select[rawCResults, !TrueQ[#["OK"]] &];

Print["Raw complex capacitance overlay: parsed OK = ", Length[goodRawC], " / ", Length[rawCResults]];

If[Length[badRawC] > 0,
	Print["Failures:"];
	Scan[
		(Print["\t", FileNameTake[#["File"]], " : ", #["Message"]]) &,
		badRawC
	];
];

If[Length[goodRawC] == 0,
	Print["No raw complex capacitance spectra could be parsed."];
	Abort[];
];

rawCMagTraces   = goodRawC[[All, "MagTrace"]];
rawCPhaseTraces = goodRawC[[All, "PhaseTrace"]];
rawCReTraces    = goodRawC[[All, "ReTrace"]];
rawCImTraces    = goodRawC[[All, "ImTrace"]];
rawCLabels      = goodRawC[[All, "Label"]];
rawCColors      = ColorData["Rainbow"] /@ Rescale[Range[Length[goodRawC]]];

(* ---------- Plots ---------- *)

rawCMagPlot =
	ListLogLogPlot[
		rawCMagTraces,
		PlotStyle -> rawCColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "|Ccomplex| (\[Micro]F)"},
		ImageSize -> 500,
		Background -> White,
		PlotLabel -> Row[{"RAW |Ccomplex|   n=", Length[goodRawC]}]
	];

rawCPhasePlot =
	ListLogLinearPlot[
		rawCPhaseTraces,
		PlotStyle -> rawCColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "Phase(Ccomplex) (deg)"},
		ImageSize -> 500,
		Background -> White,
		PlotLabel -> Row[{"RAW Phase(Ccomplex)   n=", Length[goodRawC]}]
	];

rawCRePlot =
	ListLogLinearPlot[
		rawCReTraces,
		PlotStyle -> rawCColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "Re(Ccomplex) (\[Micro]F)"},
		ImageSize -> 500,
		Background -> White,
		PlotLabel -> "RAW Re(Ccomplex)"
	];

rawCImPlot =
	ListLogLinearPlot[
		rawCImTraces,
		PlotStyle -> rawCColors,
		Joined -> True,
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Frequency (Hz)", "Im(Ccomplex) (\[Micro]F)"},
		ImageSize -> 500,
		Background -> White,
		PlotLabel -> "RAW Im(Ccomplex)"
	];

GraphicsGrid[
	{
		{rawCMagPlot},
		{rawCPhasePlot},
		{rawCRePlot},
		{rawCImPlot}
	},
	ImageSize -> 600,
	Spacings -> {0, 0}
]


rawCImPlot


(* ============================================================ *)
(* Nyquist plots for RAW and FIT impedance                      *)
(* ============================================================ *)

(* ---------- RAW Z traces DIRECTLY from imported Z ---------- *)

rawZResults =
	Table[
		Module[
			{
				file, raw, zRaw, trace, ok
			},

			file = files[[k]];
			raw = Quiet @ Check[DCT`Private`importEISTxt[file], $Failed];

			If[
				raw === $Failed || FailureQ[raw] || !AssociationQ[raw] ||
				!KeyExistsQ[raw, "Z"],

				<|
					"OK" -> False,
					"File" -> file
				|>,

				zRaw = raw["Z"];

				trace =
					Transpose[
						{
							Re[zRaw],
							-Im[zRaw]
						}
					];

				ok =
					VectorQ[
						trace,
						MatchQ[#, {_?NumericQ, _?NumericQ}] &
					];

				<|
					"OK" -> ok,
					"File" -> file,
					"Trace" -> trace
				|>
			]
		],
		{k, Length[files]}
	];

goodRawZ = Select[rawZResults, TrueQ[#["OK"]] &];

rawZTraces = goodRawZ[[All, "Trace"]];

(* ---------- FIT Z traces ---------- *)

fitZTraces =
	Table[
		Transpose[
			{
				Re[1/goodFullFitY[[k, "YFit"]]],
				-Im[1/goodFullFitY[[k, "YFit"]]]
			}
		],
		{k, Length[goodFullFitY]}
	];

(* ---------- Shared axis range ---------- *)

allNyquistVals =
	Select[
		Join[
			Flatten[rawZTraces[[All, All, 1]]],
			Flatten[rawZTraces[[All, All, 2]]],
			Flatten[fitZTraces[[All, All, 1]]],
			Flatten[fitZTraces[[All, All, 2]]]
		],
		NumericQ
	];

nyquistMax = Max[allNyquistVals];

Print["Nyquist max = ", nyquistMax];

scaleFactor = 1;

(* ---------- RAW Nyquist ---------- *)

rawNyquistPlot =
	ListLinePlot[
		rawZTraces,
		PlotStyle -> rawYColors,
		Frame -> True,
		Axes -> False,
		AspectRatio -> 1,
		FrameLabel -> {
			"Re(Z) (\[CapitalOmega])",
			"-Im(Z) (\[CapitalOmega])"
		},
		PlotRange -> {
			{0, nyquistMax/scaleFactor},
			{0, nyquistMax/scaleFactor}
		},
		ImageSize -> 700,
		Background -> White,
		PlotLabel -> Row[{"RAW Nyquist   n=", Length[goodRawZ]}]
	];

(* ---------- FIT Nyquist ---------- *)

fitNyquistPlot =
	ListLinePlot[
		fitZTraces,
		PlotStyle -> fitYColors,
		Frame -> True,
		Axes -> False,
		AspectRatio -> 1,
		FrameLabel -> {
			"Re(Z) (\[CapitalOmega])",
			"-Im(Z) (\[CapitalOmega])"
		},
		PlotRange -> {
			{0, nyquistMax/scaleFactor},
			{0, nyquistMax/scaleFactor}
		},
		ImageSize -> 700,
		Background -> White,
		PlotLabel -> Row[{"FIT Nyquist   n=", Length[goodFullFitY]}]
	];

GraphicsRow[
	{
		rawNyquistPlot,
		fitNyquistPlot
	},
	ImageSize -> 1400
]


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ============================== *)
(* User control: area window      *)
(* ============================== *)
fMinArea = 10/(2 Pi);
fMaxArea = 500/(2 Pi);
fMid = 100/(2 Pi);   (* Hz, user-selected split point *)

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
(* ---------- helpers for moments over log10(k) ---------- *)

areaOverLog10k[k_List, g_List] := Module[{x},
  If[Length[k] < 2 || Length[g] != Length[k], Return[0.0]];
  x = Log10[k];
  Total[
    Table[
      0.5 (g[[j]] + g[[j + 1]]) (x[[j + 1]] - x[[j]]),
      {j, 1, Length[k] - 1}
    ]
  ]
];

firstMomentOverLog10k[k_List, g_List] := Module[{x, area, num},
  If[Length[k] < 2 || Length[g] != Length[k], Return[Missing["TooFewPoints"]]];
  x = Log10[k];
  area = areaOverLog10k[k, g];
  If[!NumericQ[area] || area <= 0, Return[Missing["BadArea"]]];
  num = Total[
    Table[
      0.5 (g[[j]] k[[j]] + g[[j + 1]] k[[j + 1]]) (x[[j + 1]] - x[[j]]),
      {j, 1, Length[k] - 1}
    ]
  ];
  num/area
];

secondCentralMomentOverLog10k[k_List, g_List] := Module[{x, area, mu1, num},
  If[Length[k] < 2 || Length[g] != Length[k], Return[Missing["TooFewPoints"]]];
  x = Log10[k];
  area = areaOverLog10k[k, g];
  If[!NumericQ[area] || area <= 0, Return[Missing["BadArea"]]];
  mu1 = firstMomentOverLog10k[k, g];
  If[!NumericQ[mu1], Return[Missing["BadFirstMoment"]]];
  num = Total[
    Table[
      0.5 (g[[j]] (k[[j]] - mu1)^2 + g[[j + 1]] (k[[j + 1]] - mu1)^2) (x[[j + 1]] - x[[j]]),
      {j, 1, Length[k] - 1}
    ]
  ];
  num/area
];

sigmaOverLog10k[k_List, g_List] := Module[{m2},
  m2 = secondCentralMomentOverLog10k[k, g];
  If[NumericQ[m2] && m2 >= 0, Sqrt[m2], Missing["BadSecondMoment"]]
];

(* ---------- main table ---------- *)

areaResultsFull = Table[
  Module[
    {
      spec, file, k, g, ord, kk, gg, keep, kkUse, ggUse, idx,
      kSlow, gSlow, kFast, gFast,
      cSlow, cFast, cTot, frac,
      mu1Slow, mu1Fast, mu2Slow, mu2Fast, sigmaSlow, sigmaFast
    },

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

    kSlow = If[idx >= 2, kkUse[[;; idx]], {}];
    gSlow = If[idx >= 2, ggUse[[;; idx]], {}];

    kFast = If[idx + 1 <= Length[kkUse] - 1, kkUse[[idx + 1 ;;]], {}];
    gFast = If[idx + 1 <= Length[kkUse] - 1, ggUse[[idx + 1 ;;]], {}];

    cSlow = If[Length[kSlow] >= 2, areaOverLog10k[kSlow, gSlow], 0.0];
    cFast = If[Length[kFast] >= 2, areaOverLog10k[kFast, gFast], 0.0];
    cTot = cSlow + cFast;
    frac = If[cTot > 0, cFast/cTot, Indeterminate];

    mu1Slow = If[Length[kSlow] >= 2, firstMomentOverLog10k[kSlow, gSlow], Missing["TooFewPoints"]];
    mu1Fast = If[Length[kFast] >= 2, firstMomentOverLog10k[kFast, gFast], Missing["TooFewPoints"]];

    mu2Slow = If[Length[kSlow] >= 2, secondCentralMomentOverLog10k[kSlow, gSlow], Missing["TooFewPoints"]];
    mu2Fast = If[Length[kFast] >= 2, secondCentralMomentOverLog10k[kFast, gFast], Missing["TooFewPoints"]];

    sigmaSlow = If[Length[kSlow] >= 2, sigmaOverLog10k[kSlow, gSlow], Missing["TooFewPoints"]];
    sigmaFast = If[Length[kFast] >= 2, sigmaOverLog10k[kFast, gFast], Missing["TooFewPoints"]];

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

      "kMeanSlow" -> mu1Slow,
      "kMeanFast" -> mu1Fast,

      "kVarSlow" -> mu2Slow,
      "kVarFast" -> mu2Fast,

      "kSigmaSlow" -> sigmaSlow,
      "kSigmaFast" -> sigmaFast,

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
      PlotRange -> {{0, cMax}, {0, 1}}
    ],
    Plot[
      langmuirFB[c, KD],
      {c, 0,  cMax},
      PlotStyle -> {GrayLevel[0.35], Thick}
    ],
    Plot[
      0.58 langmuirFB[c, 70] + 0.08,
      {c, 0,  cMax},
      PlotStyle -> {GrayLevel[0.35], Dashed}
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

fbTimeData = ({#["TimeHr"], #["FractionBound"]} & /@ fbTimeRows)

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


fbTimeRows = SortBy[
  Select[
    areaResultsFull,
    NumericQ[#["TimeHr"]] &&
    NumericQ[#["cSlow_F"]] &&
    NumericQ[#["cFast_F"]] &&
    NumericQ[#["kMeanSlow"]] &&
    NumericQ[#["kMeanFast"]] &&
    NumericQ[#["kVarSlow"]] &&
    NumericQ[#["kVarFast"]] &
  ],
  #["TimeHr"] &
];

fbSlowTimeData = ({#["TimeHr"], #["cSlow_F"]} & /@ fbTimeRows);
fbFastTimeData = ({#["TimeHr"], #["cFast_F"]} & /@ fbTimeRows);

fbSlowMeanData = ({#["TimeHr"], #["kMeanSlow"]} & /@ fbTimeRows);
fbFastMeanData = ({#["TimeHr"], #["kMeanFast"]} & /@ fbTimeRows);

fbSlowVarData = ({#["TimeHr"], #["kVarSlow"]} & /@ fbTimeRows);
fbFastVarData = ({#["TimeHr"], #["kVarFast"]} & /@ fbTimeRows);

slowStyle = Directive[RGBColor[0.8, 0.2, 0.2], Thick];
fastStyle = Directive[RGBColor[0.15, 0.6, 0.35], Thick];

capacitancePlot =
 ListPlot[
  {fbSlowTimeData, fbFastTimeData},

  Joined -> True,
  PlotStyle -> {slowStyle, fastStyle},
  PlotMarkers -> None,

  PlotLegends -> Placed[
    LineLegend[
      {slowStyle, fastStyle},
      {"C_slow", "C_fast"}
    ],
    Right
  ],

  Background -> White,
  PlotRangePadding -> Scaled[0.02],

  Frame -> True,
  FrameStyle -> Directive[Black, AbsoluteThickness[1.2]],

  FrameLabel -> {
    Style["Time (hours)", 18],
    Style["Capacitance (F)", 18]
  },
  BaseStyle -> {FontFamily -> "Arial", 14},
  LabelStyle -> Directive[18],
  FrameTicksStyle -> Directive[14],
  FrameTicks -> {
    {Automatic, None},
    {Automatic, None}
  },
  PlotRange -> {All, {0, 2 10^-7}},
  ImageSize -> 600,
  PlotLabel -> Style["Fast and slow state capacitance vs time", 16]
 ];

firstMomentPlot =
 ListPlot[
  {fbSlowMeanData, fbFastMeanData},

  Joined -> True,
  PlotStyle -> {slowStyle, fastStyle},
  PlotMarkers -> None,

  PlotLegends -> Placed[
    LineLegend[
      {slowStyle, fastStyle},
      {"k_mean, slow", "k_mean, fast"}
    ],
    Right
  ],

  Background -> White,
  PlotRangePadding -> Scaled[0.02],

  Frame -> True,
  FrameStyle -> Directive[Black, AbsoluteThickness[1.2]],

  FrameLabel -> {
    Style["Time (hours)", 18],
	Style[Row[{"First moment of k (", Superscript["s", -1], ")"}], 18]
  },

  BaseStyle -> {FontFamily -> "Arial", 14},
  LabelStyle -> Directive[18],
  FrameTicksStyle -> Directive[14],

  FrameTicks -> {
    {Automatic, None},
    {Automatic, None}
  },

  PlotRange -> {Automatic,{0,350}},
  ImageSize -> 600,

  PlotLabel -> Style["Fast and slow state first moment vs time", 16]
 ];

secondMomentPlot =
 ListPlot[
  {fbSlowVarData, fbFastVarData},

  Joined -> True,
  PlotStyle -> {slowStyle, fastStyle},
  PlotMarkers -> None,

  PlotLegends -> Placed[
    LineLegend[
      {slowStyle, fastStyle},
      {"k_var, slow", "k_var, fast"}
    ],
    Right
  ],

  Background -> White,
  PlotRangePadding -> Scaled[0.02],

  Frame -> True,
  FrameStyle -> Directive[Black, AbsoluteThickness[1.2]],

  FrameLabel -> {
    Style["Time (hours)", 18],
    Style["Second central moment of k ((s^-1)^2)", 18]
  },

  BaseStyle -> {FontFamily -> "Arial", 14},
  LabelStyle -> Directive[18],
  FrameTicksStyle -> Directive[14],

  FrameTicks -> {
    {Automatic, None},
    {Automatic, None}
  },

  PlotRange -> {Automatic,{0,500}},
  ImageSize -> 600,

  PlotLabel -> Style["Fast and slow state second moment vs time", 16]
 ];

Column[
  {
    capacitancePlot,
    firstMomentPlot,
    secondMomentPlot
  },
  Spacings -> 2
]


(* ============================================================ *)
(* TIME-SERIES / PSD ANALYSIS OF EIS FOURIER COMPONENTS         *)
(* Re-imports raw Z directly from files in goodSpecs             *)
(* Last 20 hours only                                           *)
(* ============================================================ *)

If[Length[goodSpecs] == 0,
	Print["goodSpecs is empty."];
	Abort[];
];

(* ---------- helper to pull finish times ---------- *)

finishTimesS = Lookup[goodSpecs[[All, "Spec"]], "FinishTimeS", Missing["NoTime"]];

validTimeIdx =
	Select[
		Range[Length[goodSpecs]],
		NumericQ[finishTimesS[[#]]] &
	];

If[Length[validTimeIdx] < 20,
	Print["Too few spectra with valid FinishTimeS."];
	Abort[];
];

timesSAll = finishTimesS[[validTimeIdx]];
specSubsetAll = goodSpecs[[validTimeIdx]];

tMax = Max[timesSAll];
tMinKeep = tMax - 20.*3600.;

keepLast20h = Map[# >= tMinKeep &, timesSAll];

timesS = Pick[timesSAll, keepLast20h];
specSubset = Pick[specSubsetAll, keepLast20h];

If[Length[specSubset] < 20,
	Print["Too few spectra in last 20 hours."];
	Abort[];
];

(* ---------- sort by time ---------- *)

ord = Ordering[timesS];
timesS = timesS[[ord]];
specSubset = specSubset[[ord]];

timesHr = (timesS - Min[timesS])/3600.0;
nSpec = Length[specSubset];

Print["Using ", nSpec, " spectra from last 20 hours."];

(* ---------- re-import raw files directly ---------- *)

rawImportedSubset =
	Table[
		DCT`Private`importEISTxt[specSubset[[k, "File"]]],
		{k, 1, nSpec}
	];

If[!VectorQ[rawImportedSubset, AssociationQ],
	Print["At least one raw file failed to import."];
	Abort[];
];

(* ---------- common frequency grid check ---------- *)

freqTemplate = rawImportedSubset[[1, "FreqHz"]];

sameFreqGridQ =
	And @@ Table[
		rawImportedSubset[[k, "FreqHz"]] === freqTemplate,
		{k, 2, nSpec}
	];

If[!sameFreqGridQ,
	Print["Frequency grids are not identical across spectra."];
	Abort[];
];

freqList = N[freqTemplate];
nFreq = Length[freqList];

(* ---------- build matrices ---------- *)

zMat =
	Table[
		rawImportedSubset[[k, "Z"]],
		{k, 1, nSpec}
	];

If[!MatrixQ[zMat, NumericQ],
	Print["zMat is not numeric. Check imported Z values."];
	Abort[];
];

rsVec =
	Table[
		specSubset[[k, "Spec", "Rs"]],
		{k, 1, nSpec}
	];

If[!VectorQ[rsVec, NumericQ],
	Print["Could not extract numeric Rs values from goodSpecs."];
	Abort[];
];

yMat =
	Table[
		1/(zMat[[k]] - rsVec[[k]]),
		{k, 1, nSpec}
	];

If[!MatrixQ[yMat, NumericQ],
	Print["yMat is not numeric. Check zMat and Rs values."];
	Abort[];
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
	{x0, n, fs, fft, kmax, psd, ff},

	n = Length[x];

	If[n < 4,
		Return[<|"PSDfreq" -> {}, "PSD" -> {}|>]
	];

	If[!VectorQ[x, NumericQ],
		Return[<|"PSDfreq" -> {}, "PSD" -> {}|>]
	];

	fs = 1./dt;
	x0 = N[x - Mean[x]];

	fft = Fourier[x0, FourierParameters -> {1, -1}];

	kmax = Floor[n/2];

	ff = Range[0, kmax] * fs/n;
	psd = (1./(fs*n)) * Abs[fft[[1 ;; kmax + 1]]]^2;

	If[kmax >= 2,
		psd[[2 ;; -2]] = 2.0 * psd[[2 ;; -2]];
	];

	<|"PSDfreq" -> ff, "PSD" -> psd|>
];

(* ---------- compute PSD and statistics at each EIS frequency ---------- *)

psdResults =
	Table[
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

			If[!VectorQ[zSeries, NumericQ] || !VectorQ[ySeries, NumericQ],
				Return[<|"OK" -> False, "Freq" -> f|>]
			];

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
			varZMagPSD =
				If[
					Length[psdZMag["PSDfreq"]] > 2,
					Total[Most[psdZMag["PSD"]] * Differences[psdZMag["PSDfreq"]]],
					Missing["NoPSD"]
				];

			meanYMag = Mean[yMagSeries];
			varYMagRaw = Variance[yMagSeries];
			varYMagPSD =
				If[
					Length[psdYMag["PSDfreq"]] > 2,
					Total[Most[psdYMag["PSD"]] * Differences[psdYMag["PSDfreq"]]],
					Missing["NoPSD"]
				];

			meanZRe = Mean[zReSeries];
			meanNegZIm = Mean[zImSeries];
			varZReRaw = Variance[zReSeries];
			varNegZImRaw = Variance[zImSeries];

			meanYRe = Mean[yReSeries];
			meanNegYIm = Mean[yImSeries];
			varYReRaw = Variance[yReSeries];
			varNegYImRaw = Variance[yImSeries];

			<|
				"OK" -> True,
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

goodPSDResults =
	Select[
		psdResults,
		AssociationQ[#] &&
			TrueQ[Lookup[#, "OK", False]] &&
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

psdCurves =
	Transpose[{#["PSDfreq"]*3600, #["PSD"]}] & /@ goodPSDResults;

psdLabels =
	Table[
		ToString@NumberForm[goodPSDResults[[k, "Freq"]], {7, 2}] <> " Hz",
		{k, Length[goodPSDResults]}
	];

psdColors = ColorData["Rainbow"] /@ Rescale[Range[Length[psdCurves]]];

psdOverlayPlot =
	ListLogLogPlot[
		psdCurves,
		PlotStyle -> psdColors,
		Frame -> True,
		Joined -> True,
		Axes -> False,
		FrameLabel -> {"Temporal frequency (cycles/hour)", "PSD of |Z|"},
		PlotRange -> All,
		ImageSize -> 900,
		Background -> White,
		PlotLegends ->
			Placed[
				LineLegend[
					psdColors,
					psdLabels,
					LegendLayout -> "Column",
					LabelStyle -> Directive[11]
				],
				Right
			],
		PlotLabel -> "Overlay PSD of |Z| for each EIS frequency"
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
yShotLine = makeReferenceLine[yVarMeanRaw, 1.5];
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





 (* ------------------------------------------------------------ *)
(* Final admittance variance plot with labels and slope guides  *)
(* ------------------------------------------------------------ *)

ySlopeList = {0.5, 1.0, 1.5, 2.0};

ySlopeLines =
	makeReferenceLine[yBasePoints, #] & /@ ySlopeList;

ySlopeLabels =
	("Slope " <> ToString[NumberForm[#, {3, 1}]] <>
	 "   Var[Y] \[Proportional] |Y|^" <> ToString[NumberForm[#, {3, 1}]]) & /@
		ySlopeList;

yVarPlot =
	ListLogLogPlot[
		Join[
			{
				yBasePoints
			},
			ySlopeLines,
			{
				yLabelPoints
			}
		],
		Joined -> Join[
			{False},
			ConstantArray[True, Length[ySlopeLines]],
			{False}
		],
		PlotMarkers -> Join[
			{
				{Automatic, Medium}
			},
			ConstantArray[None, Length[ySlopeLines]],
			{
				None
			}
		],
		PlotStyle -> Join[
			{
				Directive[Black]
			},
			{
				Directive[Blue, Dashed, Thick],
				Directive[Darker[Green], Dashed, Thick],
				Directive[Purple, Dashed, Thick],
				Directive[Orange, Dashed, Thick]
			},
			{
				Directive[Black]
			}
		],
		Frame -> True,
		Axes -> False,
		FrameLabel -> {"Mean |Y|", "Variance of |Y|"},
		PlotRange -> {{0.001,Automatic},{10^-9, Automatic}},
		ImageSize -> 800,
		Background -> White,
		PlotLegends -> Placed[
			Join[
				{
					"Raw variance"
				},
				ySlopeLabels,
				{
					"Labeled representative frequencies"
				}
			],
			Right
		],
		PlotLabel -> "Admittance variance vs mean (labeled)"
	];

yVarPlot


(* get concnetration by inverting Langmuir isotherm*)


