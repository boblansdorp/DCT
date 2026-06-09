(* ::Package:: *)

(* DCTSensitivity
   Parameter sensitivity sweeps for a single EIS file.
   Each function refits the data over a range of one parameter while
   holding the others fixed, returning a list of result Associations.

   Public API:
     LambdaSweep[file, lambdas, baseOpts]
       -> list of <|"Lambda", "Spec"|> Associations

     FrequencyRangeSweep[file, {fMins, fMaxs}, baseOpts]
       -> list of <|"FMin", "FMax", "Spec"|> Associations

     WeightPowerSweep[file, wPows, baseOpts]
       -> list of <|"WeightPower", "Spec"|> Associations
*)

BeginPackage["DCTSensitivity`", {"DCT`"}]

LambdaSweep::usage =
  "LambdaSweep[file, {l1,l2,...}, opts] fits file at each lambda value. \
Returns list of <|\"Lambda\"->, \"Spec\"->|>."

FrequencyRangeSweep::usage =
  "FrequencyRangeSweep[file, {{fMin1,fMax1},...}, opts] fits file at each \
frequency range. Returns list of <|\"FMin\", \"FMax\", \"Spec\"|>."

WeightPowerSweep::usage =
  "WeightPowerSweep[file, {p1,p2,...}, opts] fits file at each WeightPower. \
Returns list of <|\"WeightPower\", \"Spec\"|>."

Begin["`Private`"]

fitOne[file_, overrideOpts_, baseOpts_] := Module[{allOpts, spec},
  allOpts = Join[overrideOpts, baseOpts];
  spec = Quiet @ Check[
    DCT`DCTSpectrum[file, Sequence @@ allOpts],
    $Failed
  ];
  spec
]

LambdaSweep[file_String, lambdas_List, baseOpts : OptionsPattern[DCT`DCTSpectrum]] :=
  Table[
    <|"Lambda" -> lam, "Spec" -> fitOne[file, {"LambdaND" -> lam}, {baseOpts}]|>,
    {lam, lambdas}
  ]

FrequencyRangeSweep[file_String, ranges_List, baseOpts : OptionsPattern[DCT`DCTSpectrum]] :=
  Table[
    <|"FMin" -> r[[1]], "FMax" -> r[[2]],
      "Spec" -> fitOne[file, {"FMinUse" -> r[[1]], "FMaxUse" -> r[[2]]}, {baseOpts}]|>,
    {r, ranges}
  ]

WeightPowerSweep[file_String, powers_List, baseOpts : OptionsPattern[DCT`DCTSpectrum]] :=
  Table[
    <|"WeightPower" -> p, "Spec" -> fitOne[file, {"WeightPower" -> p}, {baseOpts}]|>,
    {p, powers}
  ]

End[]
EndPackage[]
