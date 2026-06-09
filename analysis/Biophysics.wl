(* ::Package:: *)

(* DCTBiophysics
   Symbolic derivation and fitting of the aptamer folding model.

   States
   ------
     U  — unfolded (electrochemically slow, k_peak in range 1)
     F  — folded, no ligand bound
     B  — folded AND ligand bound
     NF — permanently non-folding (optional 4th state)

   Electrochemical observables
   ---------------------------
     [F*] = [F] + [B]   (appears as the fast/high-k peak)
     [U*] = [U]          (3-state model)
     [U*] = [U] + [NF]  (4-state model)

   Equilibria
   ----------
     KS    = [F]/[U]           structural switching constant
     KDint = [F][C]/[B]        intrinsic ligand dissociation constant
                               (KD notation used throughout for brevity)

   Public API
   ----------
     $FF3State       — symbolic ff(C, KS, KD): 3-state closed-form
     $FF4State       — symbolic ff(C, KS, KD, NF): 4-state closed-form

     FitFoldedFraction[pts]
       -> Association with keys Model3, Model4, each containing
          BestFitParameters, FitFunction (pure function of C),
          and the FitObject itself.

     PlotFoldedFraction[pts, fitResult]
       -> Graphics overlaying data and both model curves.
*)

BeginPackage["DCTBiophysics`"]

$FF3State::usage =
  "Symbolic folded fraction ff(C, KS, KD) for the 3-state model \
(U, F, B) without a non-folding population."

$FF4State::usage =
  "Symbolic folded fraction ff(C, KS, KD, NF) for the 4-state model \
(U, F, B, NF) with a non-folding population NF."

FitFoldedFraction::usage =
  "FitFoldedFraction[pts] fits both 3-state and 4-state models to \
{{concentration, ff}, ...} data. Returns an Association."

PlotFoldedFraction::usage =
  "PlotFoldedFraction[pts, fitResult] plots data and model overlays \
on a log-concentration axis."

conc::usage = "Concentration variable in $FF3State and $FF4State."
KS::usage   = "Structural switching constant in $FF3State and $FF4State."
KD::usage   = "Ligand dissociation constant in $FF3State and $FF4State."
NF::usage   = "Non-folding fraction in $FF4State."

Begin["`Private`"]

(* ================================================================ *)
(* Symbolic derivation                                               *)
(* ================================================================ *)

(* --- 3-state: U + F + B = 1 --- *)
$FF3State = Module[
  {uSym, fSym, bSym, sol, fstar, ustar},

  (* equilibrium relations *)
  (* fSym = KS * uSym  =>  fSym - KS uSym = 0 *)
  (* bSym = fSym * C / KD  =>  bSym - (KS uSym C)/KD = 0 *)
  (* normalisation *)
  sol = Solve[
    {fSym == KS * uSym,
     bSym == fSym * conc / KD,
     uSym + fSym + bSym == 1},
    {uSym, fSym, bSym}
  ][[1]];

  fstar = (fSym + bSym) /. sol;
  ustar = uSym /. sol;
  With[{e = Simplify[fstar / (fstar + ustar)]},
    Numerator[e] / Collect[Denominator[e], KS]]
];

(* --- 4-state: U + F + B + NF = 1 --- *)
$FF4State = Module[
  {uSym, fSym, bSym, sol, fstar, ustar},

  (* NF fraction is fixed; the remaining (1-NF) is distributed over U/F/B *)
  sol = Solve[
    {fSym == KS * uSym,
     bSym == fSym * conc / KD,
     uSym + fSym + bSym == 1 - NF},
    {uSym, fSym, bSym}
  ][[1]];

  fstar = (fSym + bSym) /. sol;
  ustar = (uSym + NF) /. sol;   (* NF contributes to the slow/unfolded observable *)
  (* Simplify can produce -(expr*(-1+NF)); fold the outer -1 in to get (1-NF) form *)
  With[{e = Simplify[fstar / (fstar + ustar)] //.
      Times[-1, a___, Plus[-1, x_Symbol], b___] :> Times[a, 1 - x, b]},
    Numerator[e] / Collect[Denominator[e], KS]]
];

(* ================================================================ *)
(* Fitting                                                            *)
(* ================================================================ *)

FitFoldedFraction[pts_List] := Module[
  {model3, model4, fit3, fit4, p3, p4, fn3, fn4},

  (* substitute symbolic expressions with numeric variables *)
  model3 = $FF3State /. {conc -> cx, KS -> ks, KD -> kd};
  model4 = $FF4State /. {conc -> cx, KS -> ks, KD -> kd, NF -> nf};

  fit3 = NonlinearModelFit[pts, model3,
    {{ks, 0.1, 0., Infinity}, {kd, 100., 0., Infinity}},
    cx,
    MaxIterations -> 500
  ];

  fit4 = NonlinearModelFit[pts, model4,
    {{ks, 0.1, 0., Infinity}, {kd, 100., 0., Infinity},
     {nf, 0.3, 0., 0.99}},
    cx,
    MaxIterations -> 500
  ];

  p3 = fit3["BestFitParameters"];
  p4 = fit4["BestFitParameters"];

  fn3 = With[{e = model3 /. p3}, (e /. cx -> #) &];
  fn4 = With[{e = model4 /. p4}, (e /. cx -> #) &];

  <|
    "Model3" -> <|
      "BestFitParameters" -> p3,
      "KS"  -> (ks /. p3),
      "KD"  -> (kd /. p3),
      "FitFunction" -> fn3,
      "FitObject"   -> fit3
    |>,
    "Model4" -> <|
      "BestFitParameters" -> p4,
      "KS"  -> (ks /. p4),
      "KD"  -> (kd /. p4),
      "NF"  -> (nf /. p4),
      "FitFunction" -> fn4,
      "FitObject"   -> fit4
    |>
  |>
]

(* ================================================================ *)
(* Plot                                                               *)
(* ================================================================ *)

PlotFoldedFraction[pts_List, fitResult_Association,
    opts : OptionsPattern[]] := Module[
  {darkBg, fgSt, cMin, cMax, fn3, fn4, ks3, kd3, ks4, kd4, nf4},

  darkBg = GrayLevel[0.12];
  fgSt   = Directive[White, AbsoluteThickness[1.1]];

  cMin = Max[0.1, Min[pts[[All,1]]]];
  cMax = Max[pts[[All,1]]];

  fn3  = fitResult["Model3", "FitFunction"];
  fn4  = fitResult["Model4", "FitFunction"];
  ks3  = fitResult["Model3", "KS"];
  kd3  = fitResult["Model3", "KD"];
  ks4  = fitResult["Model4", "KS"];
  kd4  = fitResult["Model4", "KD"];
  nf4  = fitResult["Model4", "NF"];

  Show[
    ListLinePlot[pts,
      ScalingFunctions -> {"Log10", None},
      Joined -> True, PlotMarkers -> {Automatic, 8},
      PlotStyle -> Directive[White, Opacity[0.7], AbsoluteThickness[1.5]],
      Frame -> True, Axes -> False,
      Background -> darkBg, FrameStyle -> fgSt,
      LabelStyle -> Directive[White, 14, FontFamily -> "Arial"],
      FrameLabel -> {Style["[Ligand] (\[Mu]M)", 14, White],
                     Style["Fraction folded", 14, White]},
      PlotRange -> {{cMin * 0.8, cMax * 1.2}, {0, 1}},
      PlotRangePadding -> Scaled[0.04],
      ImageSize -> 520, AspectRatio -> 0.75
    ],
    Plot[fn3[cx], {cx, cMin, cMax},
      ScalingFunctions -> {"Log10", None},
      PlotStyle -> Directive[GrayLevel[0.65], Dashed, AbsoluteThickness[2.]],
      Axes -> False],
    Plot[fn4[cx], {cx, cMin, cMax},
      ScalingFunctions -> {"Log10", None},
      PlotStyle -> Directive[White, AbsoluteThickness[2.]],
      Axes -> False],
    Epilog -> {
      Text[Style[Row[{"3-state:  K", Subscript["S",""], " = ",
        NumberForm[ks3, {3,2}], ",  K", Subscript["D",""],
        " = ", NumberForm[Round[kd3, 1.], {5,1}], " \[Mu]M"}],
        GrayLevel[0.65], 12, FontFamily -> "Arial"], Scaled[{0.55, 0.93}]],
      Text[Style[Row[{"4-state:  K", Subscript["S",""], " = ",
        NumberForm[ks4, {3,2}], ",  K", Subscript["D",""],
        " = ", NumberForm[Round[kd4, 1.], {5,1}], " \[Mu]M",
        ",  NF = ", NumberForm[nf4, {3,2}]}],
        White, 12, FontFamily -> "Arial"], Scaled[{0.50, 0.50}]]
    },
    opts
  ]
]

End[]
EndPackage[]
