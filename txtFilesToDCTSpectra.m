(* ::Package:: *)

(* ============================================================ *)
(* DCT batch runner                                              *)
(*   1) Find all .txt files in a directory                        *)
(*   2) Fit each file to a DCT spectrum                            *)
(*   3) Plot the results                                           *)
(* ============================================================ *)


(* ---------- LOAD PACKAGE ---------- *)
ClearAll["DCT`*"]; (* first unload it *)
dctPackagePath = "C:\\Users\\Bob Lansdorp\\Documents\\DCT\\DCT.wl";
Get[dctPackagePath]

ClearAll["NNLS`*"]; (* first unload it *)
NNLSPackagePath = "C:\\Users\\Bob Lansdorp\\Documents\\DCT\\NNLSFit.m";
Get[NNLSPackagePath]






(* ---------- USER SETTINGS ---------- *)



dataDir = "C:\\Users\\Bob Lansdorp\\Documents\\DCT\\data\\test";
dataDir = "C:\\Users\\Bob Lansdorp\\Documents\\DCT\\data\\2026-02-25";
lambdaND = 1 10^-4;

debugFlag = False;
fileDecimation = 1;   (* keep every Nth file: 10 -> ~450/10 = 45 files *)



(* fMinUse = 0.1; *)     (* Hz *)
fMinUse = 10;      (* Hz *) (* starting to see resistive behavior at low freq (oxygen reduction? diffusion?) *)
fMaxUse = 800;     (* Hz *)

binsPerDecade = 25;

weightPower = 1.0; (* how much do we weight each data point? around 0.5 or 1 works, has to do with SNR of potentiostat *)


(* Optional: restrict the frequency range (must match package options) *)
paddingDecades = 1.0; (* sets how many decades beyond the measured frequency range the tau values extend *)
(* 0 = no padding, 1 = a decdade of apadding. *)

tauMinFactor = 10^-paddingDecades;
tauMaxFactor = 10^paddingDecades;          (* tauMaxUse = tauMax * factor *) (* add factor of 10 to the maximum as padding *)



(* ============================================================ *)
(* 1) FIND ALL TXT FILES                                         *)
(* ============================================================ *)
txtFilesAll = Sort @ FileNames["*.txt", dataDir];

Print["Found ", Length[txtFilesAll], " .txt files total."];
If[Length[txtFilesAll] == 0, Abort[]];

fileDecimation = Max[1, Round[fileDecimation]];

txtFiles = txtFilesAll[[;; ;; fileDecimation]];

Print["Using ", Length[txtFiles], " files after decimation (every ", fileDecimation, "th file)."];



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
            "WeightPower" -> weightPower
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
    PlotRange -> {Automatic, {0, 5 10^-6}},
    PlotLegends -> Placed[labels, Right],
    ImageSize -> 700
  ]
]


(* ::InheritFromParent:: *)
(**)


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
  Frame -> True,

  FrameLabel -> {
	Style["Time (hours)", lblSize],
    Style["k (s^-1)", lblSize]
  },

  (* Log axis on k *)
  ScalingFunctions -> {None, "Log10"},

  (* Make ticks readable *)
  BaseStyle -> {FontFamily -> "Arial", tickSize},
  FrameTicksStyle -> Directive[tickSize],
  LabelStyle -> Directive[lblSize],

  (* Less clutter: fewer ticks on x, decent ticks on log y *)
  FrameTicks -> {
  {Automatic, None},
  {Automatic, None}
},

  PlotLegends -> Placed[Automatic, Right],

  PlotRange -> All,
  InterpolationOrder -> 0,     (* pixel look; change to 1 for smoother *)
  ImageSize -> 900,

  PlotLabel -> Style[
    Row[{
      "g(k) heat map   nExp=", nExp,
      "   k-range=[", ScientificForm[kMin, 3], ", ", ScientificForm[kMax, 3], "] s^-1"
    }],
    titleSize
  ]
]



(* ::InheritFromParent:: *)
(**)


(* ============================================================ *)
(* Import raw EIS txt (Gamry-style semicolon separated)          *)
(* Works for:                                                    *)
(*   Index;Frequency (Hz);Z'...                                  *)
(*   Frequency (Hz);Z'...                                        *)
(* Returns <|"FreqHz"->..., "Z"->...|>                            *)
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
(* Nyquist sanity plot: raw-from-file points + spec fit line     *)
(* ============================================================ *)

selector = 41;

selector = Clip[selector, {1, Length[goodSpecs]}];
file = goodSpecs[[selector, "File"]];
spec = goodSpecs[[selector, "Spec"]];

raw = ImportEISRaw[file];
If[raw === $Failed,
  Print["Could not parse raw file: ", file];
  Abort[];
];

fRaw = raw["FreqHz"];
zRaw = raw["Z"];

fUse = spec["FreqHz"];
zFit = spec["ZFit"];

(* Print both ranges so it\[CloseCurlyQuote]s obvious what\[CloseCurlyQuote]s being compared *)
Print["Nyquist selector = ", selector, " / ", Length[goodSpecs], " : ", FileNameTake[file]];
Print["Raw f range (Hz)  = {", Min[fRaw], ", ", Max[fRaw], "}  (n=", Length[fRaw], ")"];
Print["Used f range (Hz) = {", Min[fUse], ", ", Max[fUse], "}  (n=", Length[fUse], ")"];
Print["Rs (fit, \[CapitalOmega]) = ", spec["Rs"]];

Show[
  ListPlot[
    Transpose[{Re[zRaw], -Im[zRaw]}],
    PlotMarkers -> {Automatic, 7},
    PlotStyle -> GrayLevel[0.35]
  ],
  ListLinePlot[
    Transpose[{Re[zFit], -Im[zFit]}],
    PlotStyle -> {Red, Thick}
  ],
  Frame -> True,
  FrameLabel -> {"Z' (\[CapitalOmega])", "-Z'' (\[CapitalOmega])"},
  PlotRange -> {{0,500000},{0,500000}},
  ImageSize -> 650,
  AspectRatio -> 1,
  PlotLabel -> Row[{
     FileNameTake[file], "    Rs=", NumberForm[spec["Rs"], {8, 3}],
     " \[CapitalOmega]    used f=[", Min[fUse], ", ", Max[fUse], "] Hz"
  }]
]


(* ============================================================ *)
(* Rainbow Nyquist overlay: all fitted spectra on one plot       *)
(* Optionally overlay raw points for one selected file           *)
(* ============================================================ *)

(* ---------- USER CONTROLS ---------- *)
selector = 1;                 (* which file to show raw points for, if enabled *)
showRawPoints = False;        (* True to overlay raw points for selector file *)
useOnlyFitRange = True;       (* True: plot only ZFit from each spec *)

(* ---------- SAFETY ---------- *)
If[Length[goodSpecs] == 0,
  Print["No goodSpecs to plot."];
  Abort[];
];

selector = Clip[selector, {1, Length[goodSpecs]}];

(* ---------- Collect fit Nyquist traces ---------- *)
labels = FileNameTake /@ (goodSpecs[[All, "File"]]);

nyqFitTraces =
  Table[
    Module[{spec = goodSpecs[[k, "Spec"]], zFit},
      zFit = spec["ZFit"];
      Transpose[{Re[zFit], -Im[zFit]}]
    ],
    {k, Length[goodSpecs]}
  ];

(* ---------- Stats for annotation ---------- *)
rsList = goodSpecs[[All, "Spec", "Rs"]];
fMinList = Min /@ (goodSpecs[[All, "Spec", "FreqHz"]]);
fMaxList = Max /@ (goodSpecs[[All, "Spec", "FreqHz"]]);

Print["Nyquist overlay: nCurves = ", Length[goodSpecs]];
Print["Used f range across curves (Hz): min = ", Min[fMinList], "   max = ", Max[fMaxList]];
Print["Rs across curves (\[CapitalOmega]): min = ", Min[rsList], "   max = ", Max[rsList]];

(* ---------- Optional raw points for one file ---------- *)
rawTrace = {};
rawLabel = "";

If[showRawPoints,
  fileSel = goodSpecs[[selector, "File"]];
  raw = ImportEISRaw[fileSel];
  If[raw === $Failed,
    Print["Could not parse raw file: ", fileSel];
    rawTrace = {};
    rawLabel = " (raw parse failed)";
    ,
    zRaw = raw["Z"];
    rawTrace = Transpose[{Re[zRaw], -Im[zRaw]}];
    rawLabel = " (raw points shown for " <> FileNameTake[fileSel] <> ")";
  ];
];

(* ---------- Rainbow colors ---------- *)
colors = ColorData["Rainbow"] /@ Rescale[Range[Length[nyqFitTraces]]];

(* ---------- Plot ---------- *)
dataLogPlot = Show[
  {
    If[showRawPoints && rawTrace =!= {},
      ListLogLogPlot[
        rawTrace,
        Joined->True,
        PlotStyle -> GrayLevel[0.65],
        PlotMarkers -> {Automatic, 6}
      ],
      Nothing
    ],

    ListLogLogPlot[
      nyqFitTraces,
      PlotStyle -> colors,
      Joined->True
    ]
  },
  Frame -> True,
  FrameLabel -> {"Z' (\[CapitalOmega])", "-Z'' (\[CapitalOmega])"},
  PlotRange -> All,
  ImageSize -> 750,
  AspectRatio -> 1,
  PlotLabel -> Row[{
     "Nyquist fit overlay (rainbow)  n=", Length[goodSpecs],
     "   Rs=[", NumberForm[Min[rsList], {8, 3}], ", ", NumberForm[Max[rsList], {8, 3}], "] \[CapitalOmega]",
     "   used f overall=[", Min[fMinList], ", ", Max[fMaxList], "] Hz",
     rawLabel
  }]
]
dataLinPlot =Show[
  {
    If[showRawPoints && rawTrace =!= {},
      ListPlot[
        rawTrace,
        Joined->True,
        PlotStyle -> GrayLevel[0.65],
        PlotMarkers -> {Automatic, 6}
      ],
      Nothing
    ],

    ListPlot[
      nyqFitTraces,
      PlotStyle -> colors,
      Joined->True
    ]
  },
  Frame -> True,
  FrameLabel -> {"Z' (\[CapitalOmega])", "-Z'' (\[CapitalOmega])"},
  PlotRange -> All,
  ImageSize -> 750,
  AspectRatio -> 1,
  PlotLabel -> Row[{
     "Nyquist fit overlay (rainbow)  n=", Length[goodSpecs],
     "   Rs=[", NumberForm[Min[rsList], {8, 3}], ", ", NumberForm[Max[rsList], {8, 3}], "] \[CapitalOmega]",
     "   used f overall=[", Min[fMinList], ", ", Max[fMaxList], "] Hz",
     rawLabel
  }]
]


(* ============================================================ *)
(* Rainbow Nyquist overlay: FULL reconstructed fit from ladder   *)
(*   Zfit(\[Omega])=Rs + 1 / (i\[Omega] C0 + Sum_k gk \[CapitalDelta]log10 (i\[Omega])/(1+i\[Omega]\[Tau]k))     *)
(* ============================================================ *)

fullFitNyqTraces =
  Table[
    Module[{spec, tau, g, c0, rs, f, \[Omega], \[CapitalDelta]log, Kmat, YintFit, Zfit},
      spec = goodSpecs[[k, "Spec"]];
      tau = spec["Tau"];
      g   = spec["g"];
      c0  = spec["C0"];
      rs  = spec["Rs"];
      f   = spec["FreqHz"];
      \[Omega]   = 2 Pi f;

      \[CapitalDelta]log = Mean[Differences[Log10[tau]]];
      Kmat = Table[(I*\[Omega][[j]])/(1 + I*\[Omega][[j]]*tau[[m]]), {j, Length[\[Omega]]}, {m, Length[tau]}];

      YintFit = (I*\[Omega])*c0 + Kmat . (g * \[CapitalDelta]log);
      Zfit    = rs + 1/YintFit;

      Transpose[{Re[Zfit], -Im[Zfit]}]
    ],
    {k, Length[goodSpecs]}
  ];

colors = ColorData["Rainbow"] /@ Rescale[Range[Length[fullFitNyqTraces]]];

fitLinPlot = Show[
  ListLinePlot[
    fullFitNyqTraces,
    PlotStyle -> colors
  ],
  Frame -> True,
  FrameLabel -> {"Z' (\[CapitalOmega])", "-Z'' (\[CapitalOmega])"},
  PlotRange -> All,
  ImageSize -> 750,
  AspectRatio -> 1,
  PlotLabel -> Row[{"FULL Maxwell ladder fit Nyquist overlay   n=", Length[goodSpecs]}]
]


fitLogPlot = Show[
  ListLogLogPlot[
    fullFitNyqTraces,
    PlotStyle -> colors,
    Joined->True
  ],
  Frame -> True,
  FrameLabel -> {"Z' (\[CapitalOmega])", "-Z'' (\[CapitalOmega])"},
  PlotRange -> All,
  Joined->True,
  ImageSize -> 750,
  AspectRatio -> 1,
  PlotLabel -> Row[{"FULL Maxwell ladder fit Nyquist overlay   n=", Length[goodSpecs]}]
]


GraphicsGrid[{{dataLinPlot, fitLinPlot},{dataLogPlot, fitLogPlot}}]
