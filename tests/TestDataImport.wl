(* ::Package:: *)

(* Tests for DCTDataImport`ImportEIS
   Run this cell in a wolfbook notebook or via:
     Get["path/to/tests/TestDataImport.wl"]
   Prints a pass/fail summary.
*)

Needs["DCTDataImport`"]

Begin["DCTTests`DataImport`Private`"]

(* ---- Minimal test harness ---- *)
$pass = 0; $fail = 0;

check[name_String, got_, expected_] := Module[{ok},
  ok = If[NumberQ[expected] && NumberQ[got],
    Abs[got - expected] / Max[Abs[expected], 10^-30] < 10^-6,
    got === expected
  ];
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

(* ---- Synthetic EIS file helper ---- *)
makeTempEIS[rows_List, header_String : ""] := Module[{file, content},
  file = CreateTemporary[];
  content = If[header =!= "",
    StringJoin[header, "\n", StringRiffle[rows, "\n"]],
    StringRiffle[rows, "\n"]
  ];
  Export[file, content, "Text"];
  file
]

(* ---- Tests ---- *)

(* T1: basic semicolon-delimited with header *)
t1File = makeTempEIS[
  {"10.0; 500.0; 100.0", "1.0; 600.0; 80.0", "0.1; 700.0; 50.0"},
  "Frequency (Hz); Z' (\[CapitalOmega]); -Z'' (\[CapitalOmega])"
];
t1 = DCTDataImport`ImportEIS[t1File];
Print @ checkTrue["T1 returns Association",      AssociationQ[t1]];
Print @ checkTrue["T1 has FreqHz key",           KeyExistsQ[t1, "FreqHz"]];
Print @ checkTrue["T1 has Z key",                KeyExistsQ[t1, "Z"]];
Print @ check["T1 row count",                    Length[t1["FreqHz"]], 3];
Print @ check["T1 sorted ascending freq[1]",     t1["FreqHz"][[1]], 0.1];
Print @ check["T1 Re[Z] at lowest freq",         Re[t1["Z"][[1]]], 700.0];
Print @ check["T1 Im[Z] negative convention",    Im[t1["Z"][[1]]], -50.0];

(* T2: tab-delimited, no header, columns {freq, Zre, -Zim} *)
t2File = makeTempEIS[
  {"5.0\t200.0\t30.0", "50.0\t180.0\t40.0"},
  ""
];
t2 = DCTDataImport`ImportEIS[t2File];
Print @ checkTrue["T2 no-header parse OK",       AssociationQ[t2]];
Print @ check["T2 row count",                    Length[t2["FreqHz"]], 2];
Print @ check["T2 Re[Z]",                        Re[t2["Z"][[1]]], 200.0];

(* T3: empty file returns Failure *)
t3File = CreateTemporary[];
Export[t3File, "", "Text"];
t3 = DCTDataImport`ImportEIS[t3File];
Print @ checkTrue["T3 empty file -> Failure",    FailureQ[t3]];

(* T4: European decimal comma is handled *)
t4File = makeTempEIS[
  {"1,0; 500,0; 100,0", "10,0; 490,0; 95,0"},
  "Frequency (Hz); Z' (\[CapitalOmega]); -Z'' (\[CapitalOmega])"
];
t4 = DCTDataImport`ImportEIS[t4File];
Print @ checkTrue["T4 decimal comma -> Association", AssociationQ[t4]];
Print @ check["T4 freq[1]",                     t4["FreqHz"][[1]], 1.0];

(* ---- Cleanup and summary ---- *)
DeleteFile /@ {t1File, t2File, t3File, t4File};

Print[""];
Print[Style["DataImport: " <> ToString[$pass] <> " passed, " <> ToString[$fail] <> " failed.",
  If[$fail == 0, Darker[Green], Darker[Red]], 14, Bold]];

End[]
