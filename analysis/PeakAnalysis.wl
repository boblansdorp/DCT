(* ::Package:: *)

(* DCTPeakAnalysis
   Peak-picking utilities for DCT spectra.
   All functions are pure: input is a spec Association or list of specs.

   Public API:
     FindDCTPeak[spec]
       -> Association <|"kPeak", "gPeak", "tauPeak", "PeakIndex"|>
          kPeak via quadratic refinement in log10(k) space.
          Returns Missing["NoPeak"] if g is all zero or flat.

     IntegrateSpectrum[spec]
       -> total distributed capacitance [F], integral of g(tau) d(log10 tau)

     PeakTable[goodSpecs]
       -> Dataset of {file, time_hr, kPeak, gPeak, cDist} per spectrum
*)

BeginPackage["DCTPeakAnalysis`"]

FindDCTPeak::usage =
  "FindDCTPeak[spec] locates the dominant peak in g(k) using quadratic \
refinement in log10(k) space. Returns <|\"kPeak\", \"gPeak\", \
\"tauPeak\", \"PeakIndex\"|> or Missing[\"NoPeak\"]."

IntegrateSpectrum::usage =
  "IntegrateSpectrum[spec] integrates g(tau) over log10(tau) to give the \
total distributed capacitance [F]."

PeakTable::usage =
  "PeakTable[goodSpecs] returns a Dataset with one row per spectrum: \
File, TimeHr, kPeak (s^-1), gPeak (F/decade), CDistributed (F)."

Begin["`Private`"]

(* ------------------------------------------------------------------ *)
(* Peak finding with quadratic refinement in log10(k) space           *)
(* ------------------------------------------------------------------ *)

FindDCTPeak[spec_Association] := Module[
  {tau, g, k, logK, iMax, left, right, x1, x2, x3, y1, y2, y3,
   a, b, logKpeak, kPeak, gPeak},

  tau = spec["Tau"];
  g   = spec["g"];

  If[!VectorQ[g, NumericQ] || Max[g] <= 0.,
    Return[Missing["NoPeak"]]
  ];

  k    = 1. / tau;
  logK = Log10[k];

  iMax = First[Ordering[g, -1]];

  (* boundary check: need at least one neighbour on each side *)
  If[iMax == 1 || iMax == Length[g],
    Return[<|
      "kPeak"    -> k[[iMax]],
      "gPeak"    -> g[[iMax]],
      "tauPeak"  -> tau[[iMax]],
      "PeakIndex" -> iMax
    |>]
  ];

  (* quadratic fit through (logK[i-1], g[i-1]), (logK[i], g[i]), (logK[i+1], g[i+1]) *)
  left  = iMax - 1;
  right = iMax + 1;
  {x1, x2, x3} = logK[[{left, iMax, right}]];
  {y1, y2, y3} = g[[{left, iMax, right}]];

  (* solve for vertex of parabola through the three points *)
  a = ((y3 - y1)/(x3 - x1) - (y2 - y1)/(x2 - x1)) / (x3 - x2);
  b = (y2 - y1) / (x2 - x1) - a * (x2 + x1);

  logKpeak = If[NumericQ[a] && a < 0,
    -b / (2 a),
    logK[[iMax]]
  ];

  kPeak = 10.^logKpeak;
  gPeak = g[[iMax]];

  <|
    "kPeak"    -> kPeak,
    "gPeak"    -> gPeak,
    "tauPeak"  -> 1. / kPeak,
    "PeakIndex" -> iMax
  |>
]

(* ------------------------------------------------------------------ *)
(* Integrate g(tau) d(log10 tau)  [trapezoidal in log space]          *)
(* ------------------------------------------------------------------ *)

IntegrateSpectrum[spec_Association] := Module[
  {tau, g, logTau, n},

  tau  = spec["Tau"];
  g    = spec["g"];

  If[!VectorQ[g, NumericQ] || Length[g] < 2, Return[0.]];

  logTau = Log10[tau];
  n = Length[g];

  Total @ Table[
    0.5 (g[[j]] + g[[j+1]]) * (logTau[[j+1]] - logTau[[j]]),
    {j, 1, n - 1}
  ]
]

(* ------------------------------------------------------------------ *)
(* Summary table for a batch of spectra                               *)
(* ------------------------------------------------------------------ *)

PeakTable[goodSpecs_List] := Module[
  {rows, headers, dataRows, fmt},

  fmt[x_?NumericQ, n_] := NumberForm[x, {n, n - 1}, ExponentFunction -> (If[Abs[#] < 3, Null, #] &)];
  fmt[other_, _]        := other;

  rows = Table[
    Module[{spec, file, peakRes, cDist, timeHr},
      spec    = s["Spec"];
      file    = s["File"];
      peakRes = FindDCTPeak[spec];
      cDist   = IntegrateSpectrum[spec];
      timeHr  = If[NumericQ[spec["FinishTimeS"]], spec["FinishTimeS"] / 3600., Missing[]];
      {
        FileNameTake[file],
        fmt[timeHr, 3],
        fmt[If[AssociationQ[peakRes], peakRes["kPeak"], Missing[]], 3],
        fmt[If[AssociationQ[peakRes], peakRes["gPeak"], Missing[]], 3],
        fmt[cDist, 3],
        fmt[spec["C0"], 3],
        fmt[spec["Rs"], 3]
      }
    ],
    {s, goodSpecs}
  ];

  headers = Style[#, Bold, White, 12] & /@ {
    "File", "Time (hr)",
    Row[{"k", Subscript["peak", ""], " (", Superscript["s", -1], ")"}],
    Row[{"g", Subscript["peak", ""], " (F/dec)"}],
    Row[{"C", Subscript["dist", ""], " (F)"}],
    Row[{"C", Subscript["0", ""], " (F)"}],
    "Rs (\[CapitalOmega])"
  };

  dataRows = Map[Style[#, White, 11] &, rows, {2}];

  Grid[
    Prepend[dataRows, headers],
    Frame          -> All,
    FrameStyle     -> GrayLevel[0.35],
    Background     -> {{}, {GrayLevel[0.22], {GrayLevel[0.12]}}},
    Dividers       -> {None, {2 -> Directive[White, AbsoluteThickness[1.]]}},
    Alignment      -> {{Left, Right, Right, Right, Right, Right, Right}, Automatic},
    Spacings       -> {1.8, 0.7},
    ItemSize       -> {{Automatic, 5, 6, 6, 6, 6, 5}, Automatic}
  ]
]

End[]
EndPackage[]
