(* ::Package:: *)

(* ============================================================ *)
(* DCT batch runner                                              *)
(*   1) Find all .txt files in a directory                        *)
(*   2) Fit each file to a DCT spectrum                            *)
(*   3) Plot the results                                           *)
(* ============================================================ *)

(* ---------- USER SETTINGS ---------- *)

dctPackagePath = "C:\\Users\\Bob Lansdorp\\Documents\\DCT\\DCT.wl";

dataDir = "C:\\Users\\Bob Lansdorp\\Documents\\DCT\\data\\test";
dataDir = "C:\\Users\\Bob Lansdorp\\Documents\\DCT\\data\\2026-02-25";
lambdaND = 1 10^-3;
debugFlag = False;


(* Optional: restrict the frequency range (must match package options) *)

(* fMinUse = 0.1; *)     (* Hz *)
fMinUse = 10;      (* Hz *) (* starting to see resistive behavior at low freq (oxygen reduction? diffusion?) *)
fMaxUse = 1000;     (* Hz *)

tauMinFactor= 2 \[Pi] 2;         (* tauMinUse = tauMin * factor *) (* since we divide frequency by 2 Pi to get tau, let's multiple back and multiply by an extra factor of four - arbitrary fudge factor! *)
tauMaxFactor = 2 \[Pi] 6;          (* tauMaxUse = tauMax * factor *)

binsPerDecade = 25;

weightPower = 0.75; (* how much do we weight each data point? around 0.5 or 1 works, has to do with SNR of potentiostat *)

fileDecimation = 10;   (* keep every Nth file: 10 -> ~450/10 = 45 files *)

(* ---------- LOAD PACKAGE ---------- *)
ClearAll["DCT`*"]; (* first unload it *)
Get[dctPackagePath]

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

(* ---------- run fits ---------- *)
results =
  Table[
    Module[{file = f, spec, msg = "", ok = True},

      spec = Quiet[
        Check[
          DCT`DCTSpectrum[
            file,
            "Debug" -> debugFlag,
            "LambdaND" -> lambdaND,
            "FMinUse" -> fMinUse,
            "FMaxUse" -> fMaxUse,
		  "TauMinFactor"  -> tauMinFactor,         (* tauMinUse = tauMin * factor *)
		  "TauMaxFactor"  -> tauMaxFactor,          (* tauMaxUse = tauMax * factor *)
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
        msg = If[spec === $Failed, "DCTSpectrum returned $Failed", "DCTSpectrum did not return a valid Association"];
      ];

      <|"File" -> file, "Spec" -> spec, "OK" -> ok, "Message" -> msg|>
    ],
    {f, txtFiles}
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

selector = 1;

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
  PlotRange -> All,
  ImageSize -> 650,
  AspectRatio -> 1,
  PlotLabel -> Row[{
     FileNameTake[file], "    Rs=", NumberForm[spec["Rs"], {8, 3}],
     " \[CapitalOmega]    used f=[", Min[fUse], ", ", Max[fUse], "] Hz"
  }]
]



