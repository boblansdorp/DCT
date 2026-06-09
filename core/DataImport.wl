(* ::Package:: *)

(* DCT`DataImport
   Parses Gamry-style EIS .txt files (semicolon, tab, or whitespace delimited).
   Detects header columns automatically; falls back to columns {1,2,3} if no header.

   Public API:
     ImportEIS[file]          -> <|"FreqHz"->{...}, "Z"->{...}, "TimeS"->{...}|>
     ImportEIS[file, True]    -> same, with debug prints
*)

BeginPackage["DCTDataImport`"]

ImportEIS::usage =
  "ImportEIS[file] imports an EIS .txt file. Returns Association with keys \
\"FreqHz\", \"Z\" (complex), \"TimeS\", sorted by ascending frequency. \
Returns a Failure object on parse error."

Begin["`Private`"]

toNumber[s_String] := Module[{ss},
  ss = StringTrim[s];
  ss = StringReplace[ss, {
    "," -> ".",
    RegularExpression["([0-9.]+)[Ee]([+-]?[0-9]+)"] -> "$1*^$2"
  }];
  Quiet @ Check[ToExpression[ss], Missing["NotNumeric"]]
]
toNumber[x_] := x

splitLine[s_String] := Module[{t = StringTrim[s]},
  Which[
    StringContainsQ[t, ";"], StringTrim /@ StringSplit[t, ";"],
    StringContainsQ[t, "\t"], StringTrim /@ StringSplit[t, "\t"],
    True, DeleteCases[StringSplit[t, Whitespace], ""]
  ]
]

detectColumn[colMap_Association, candidates_List] :=
  FirstCase[Lookup[colMap, #, Missing[]] & /@ candidates, _Integer, Missing["NoCol"]]

ImportEIS[file_String] := ImportEIS[file, False]

ImportEIS[file_String, debug_] := Module[
  {rawLines, lines, header, hasHeader,
   splitHeader, headerNorm, colMap,
   freqCol, zreCol, zimCol, timeCol,
   splitRows, numRows, good,
   freqHz, zre, zimNeg, z, timeS, ord},

  rawLines = Import[file, "Lines"];
  rawLines = Select[rawLines, StringTrim[#] =!= "" &];

  If[rawLines === {},
    Return[Failure["DCTImport", <|"Message" -> "File is empty: " <> file|>]]
  ];

  header = First[rawLines];
  hasHeader = !StringMatchQ[StringTrim[header], NumberString ~~ ___];
  lines = If[hasHeader, Rest[rawLines], rawLines];

  If[hasHeader,
    splitHeader = splitLine[header];
    headerNorm = ToLowerCase[StringReplace[
      StringTrim /@ splitHeader,
      {"\[CapitalOmega]" -> "ohm", "''" -> "doubleprime", "'" -> "prime",
       "\"" -> "", " " -> "", "(" -> "", ")" -> ""}
    ]];
    colMap = AssociationThread[headerNorm -> Range[Length[headerNorm]]];

    freqCol = detectColumn[colMap, {"frequencyhz","freq/hz","freqhz","fhz"}];
    zreCol  = detectColumn[colMap, {"zprimeohm","rez/ohm","rezohm","realzohm","zreohm"}];
    zimCol  = detectColumn[colMap, {"-zdoubleprimeohm","-imz/ohm","-imzohm","minusimzohm","-imagzohm"}];
    timeCol = detectColumn[colMap, {"times","time/s"}];

    If[debug, Print["[ImportEIS] columns: freq=", freqCol, " zre=", zreCol, " zim=", zimCol, " t=", timeCol]];

    If[MemberQ[{freqCol, zreCol, zimCol}, Missing["NoCol"]],
      Return[Failure["DCTImport",
        <|"Message" -> "Required columns not found.", "Header" -> splitHeader|>]]
    ],

    freqCol = 1; zreCol = 2; zimCol = 3; timeCol = Missing["NoCol"];
    If[debug, Print["[ImportEIS] no header, using columns 1,2,3"]]
  ];

  splitRows = splitLine /@ lines;
  numRows   = (toNumber /@ #) & /@ splitRows;
  good = Select[numRows,
    Length[#] >= Max[freqCol, zreCol, zimCol] &&
    NumericQ[#[[freqCol]]] && NumericQ[#[[zreCol]]] && NumericQ[#[[zimCol]]] &
  ];

  If[good === {},
    Return[Failure["DCTImport", <|"Message" -> "No valid numeric rows in: " <> file|>]]
  ];

  freqHz = good[[All, freqCol]];
  zre    = good[[All, zreCol]];
  zimNeg = good[[All, zimCol]];
  z      = zre - I * zimNeg;

  timeS = If[IntegerQ[timeCol] && AllTrue[good, Length[#] >= timeCol &],
    good[[All, timeCol]],
    ConstantArray[Missing["NoTime"], Length[good]]
  ];

  ord = Ordering[freqHz];
  <|"FreqHz" -> freqHz[[ord]], "Z" -> z[[ord]], "TimeS" -> timeS[[ord]]|>
]

End[]
EndPackage[]
