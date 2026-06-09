(* ::Package:: *)

(* DCTAdmittanceStats
   Computes per-frequency admittance statistics across a drift ensemble:
   mean |Y(f)|, variance |Y(f)|, coefficient of variation CV(f),
   and a power-law fit Var ~ |Y_mean|^alpha, giving weightPower = -alpha.

   Requires all files to share the same frequency grid (files that don't match
   the first file's grid are skipped with a warning).

   Public API:
     ComputeAdmittanceStats[files]
       -> Association with keys:
            FreqHz, YMean, YVar, CV,
            PowerLawAlpha, PowerLawC, WeightPowerDerived,
            FilesOK, FilesTotal

     CVValidFreqRange[stats]
     CVValidFreqRange[stats, cvThresh]
       -> {fMin, fMax} for the largest contiguous window where CV < cvThresh
*)

BeginPackage["DCTAdmittanceStats`", {"DCTDataImport`"}]

ComputeAdmittanceStats::usage =
  "ComputeAdmittanceStats[files] loads all EIS files, aligns on the first \
file's frequency grid, and returns per-frequency admittance statistics \
across the ensemble. Returns an Association or a Failure object."

CVValidFreqRange::usage =
  "CVValidFreqRange[stats] or CVValidFreqRange[stats, cvThresh] returns \
{fMin, fMax} for the largest contiguous window where CV(f) < cvThresh \
(default 0.1). Falls back to the full range if no window passes."

FilterAdmittanceStats::usage =
  "FilterAdmittanceStats[stats, fMin, fMax] restricts the stats Association \
to frequencies in [fMin, fMax] and refits the power law on that window. \
Returns an updated Association with new PowerLawAlpha, PowerLawC, and \
WeightPowerDerived."

MergeAdmittanceStats::usage =
  "MergeAdmittanceStats[{stats1, stats2, ...}] pools the per-frequency \
(YMean, YVar) data from multiple per-electrode stats Associations and refits \
the power law across all electrodes. Each electrode contributes its own \
within-electrode variance estimates; the merged fit captures the shared \
noise scaling across the full admittance range."

Begin["`Private`"]

ComputeAdmittanceStats[files_List] := Module[
  {allData, goodData, freqRef, compatible, yMagAll,
   yMean, yVar, cv, goodPairs, goodLogPts, lm, params, alpha, c},

  allData = Quiet @ Check[DCTDataImport`ImportEIS[#], $Failed] & /@ files;
  goodData = Select[allData, AssociationQ];

  If[Length[goodData] < 2,
    Return[Failure["AdmittanceStats",
      <|"Message" -> "Need at least 2 valid files, got " <> ToString[Length[goodData]]|>]]
  ];

  freqRef = goodData[[1, "FreqHz"]];

  compatible = Select[goodData,
    Length[#["FreqHz"]] == Length[freqRef] &&
    Max[Abs[#["FreqHz"] - freqRef]] < 1.*^-3 &
  ];

  If[Length[compatible] < Length[goodData],
    Message[ComputeAdmittanceStats::freqmismatch,
      Length[compatible], Length[goodData]]
  ];

  If[Length[compatible] < 2,
    Return[Failure["AdmittanceStats",
      <|"Message" -> "Fewer than 2 files share the reference frequency grid"|>]]
  ];

  (* {nFiles, nFreqs} matrix of |Y(f)| per file *)
  yMagAll = Abs[1.0 / #["Z"]] & /@ compatible;

  yMean = Mean[yMagAll];
  yVar  = Variance[yMagAll];
  cv    = Sqrt[yVar] / yMean;

  (* power-law fit in log-log space: log10(Var) = alpha*log10(|Y_mean|) + log10(c) *)
  goodPairs   = Select[Transpose[{yMean, yVar}], #[[1]] > 0 && #[[2]] > 0 &];
  goodLogPts  = {Log10[#[[1]]], Log10[#[[2]]]} & /@ goodPairs;

  {alpha, c} = If[Length[goodLogPts] >= 3,
    Module[{lmFit, p},
      lmFit = LinearModelFit[goodLogPts, {1, x}, x];
      p = lmFit["BestFitParameters"];
      {p[[2]], 10.^p[[1]]}
    ],
    {2., 1.}   (* fallback if too few clean points *)
  ];

  <|
    "FreqHz"             -> freqRef,
    "YMean"              -> yMean,
    "YVar"               -> yVar,
    "CV"                 -> cv,
    "PowerLawAlpha"      -> alpha,
    "PowerLawC"          -> c,
    "WeightPowerDerived" -> -alpha * 1.,
    "FilesOK"            -> Length[compatible],
    "FilesTotal"         -> Length[files]
  |>
]

ComputeAdmittanceStats::freqmismatch =
  "Only `1` of `2` loaded files matched the reference frequency grid; others excluded."

CVValidFreqRange[stats_Association, cvThresh_: 0.1] := Module[
  {freq, cv, good, runs, goodRuns, longest},

  freq = stats["FreqHz"];
  cv   = stats["CV"];
  good = Map[TrueQ[# < cvThresh] &, cv];

  runs     = Split[Range[Length[good]], good[[#1]] === good[[#2]] &];
  goodRuns = Select[runs, good[[First[#]]] &];

  If[goodRuns === {}, Return[{Min[freq], Max[freq]}]];

  longest = First @ MaximalBy[goodRuns, Length];
  {freq[[First[longest]]], freq[[Last[longest]]]}
]

FilterAdmittanceStats[stats_Association, fMin_?NumericQ, fMax_?NumericQ] := Module[
  {freq, mask, freqF, yMeanF, yVarF, cvF,
   goodPairs, goodLogPts, alpha, c},

  freq  = stats["FreqHz"];
  mask  = MapThread[#1 >= fMin && #1 <= fMax &, {freq}];

  freqF  = Pick[freq,           mask];
  yMeanF = Pick[stats["YMean"], mask];
  yVarF  = Pick[stats["YVar"],  mask];
  cvF    = Pick[stats["CV"],    mask];

  goodPairs  = Select[Transpose[{yMeanF, yVarF}], #[[1]] > 0 && #[[2]] > 0 &];
  goodLogPts = {Log10[#[[1]]], Log10[#[[2]]]} & /@ goodPairs;

  {alpha, c} = If[Length[goodLogPts] >= 3,
    Module[{lmFit, p},
      lmFit = LinearModelFit[goodLogPts, {1, x}, x];
      p = lmFit["BestFitParameters"];
      {p[[2]], 10.^p[[1]]}
    ],
    {2., 1.}
  ];

  Join[stats, <|
    "FreqHz"             -> freqF,
    "YMean"              -> yMeanF,
    "YVar"               -> yVarF,
    "CV"                 -> cvF,
    "PowerLawAlpha"      -> alpha,
    "PowerLawC"          -> c,
    "WeightPowerDerived" -> -alpha * 1.
  |>]
]

MergeAdmittanceStats[statsList_List] := Module[
  {good, allFreq, allYMean, allYVar, allCV,
   goodPairs, goodLogPts, alpha, c},

  good = Select[statsList, AssociationQ];
  If[Length[good] == 0,
    Return[Failure["AdmittanceStats", <|"Message" -> "No valid electrode stats to merge"|>]]
  ];

  (* Pool per-frequency data from all electrodes *)
  allFreq  = Join @@ (#["FreqHz"] & /@ good);
  allYMean = Join @@ (#["YMean"]  & /@ good);
  allYVar  = Join @@ (#["YVar"]   & /@ good);
  allCV    = Join @@ (#["CV"]     & /@ good);

  goodPairs  = Select[Transpose[{allYMean, allYVar}], #[[1]] > 0 && #[[2]] > 0 &];
  goodLogPts = {Log10[#[[1]]], Log10[#[[2]]]} & /@ goodPairs;

  {alpha, c} = If[Length[goodLogPts] >= 3,
    Module[{lmFit, p},
      lmFit = LinearModelFit[goodLogPts, {1, x}, x];
      p = lmFit["BestFitParameters"];
      {p[[2]], 10.^p[[1]]}
    ],
    {2., 1.}
  ];

  <|
    "FreqHz"             -> allFreq,
    "YMean"              -> allYMean,
    "YVar"               -> allYVar,
    "CV"                 -> allCV,
    "PowerLawAlpha"      -> alpha,
    "PowerLawC"          -> c,
    "WeightPowerDerived" -> -alpha * 1.,
    "FilesOK"            -> Total[#["FilesOK"]    & /@ good],
    "FilesTotal"         -> Total[#["FilesTotal"] & /@ good]
  |>
]

End[]
EndPackage[]
