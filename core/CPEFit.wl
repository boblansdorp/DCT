(* ::Package:: *)

(* CPEFit
   Equivalent-circuit fitting for EIS spectra (constant-phase-element model).

   Circuit:   Rs --- [ Cdl  ||  (Rct --- CPE) ]

     faradaic branch:  Z_far = Rct + 1/(Q (i w)^n)
     parallel:         Y_par = i w Cdl + 1/Z_far
     total:            Z(w)  = Rs + 1/Y_par

   Parameters: Rs (Ohm), Cdl (F), Rct (Ohm), Q (S s^n), n (dimensionless, 0<n<=1).

   The fit minimises a weighted complex least-squares objective in log-parameters
   (so Rs, Cdl, Rct, Q stay positive and span decades): FindMinimum from
   data-driven seeds, with an NMinimize / DifferentialEvolution fallback.

   Public API:
     CPEImpedance[{Rs, Cdl, Rct, Q, n}, freqHz]  -> complex Z list (forward model)
     FitCPE[freqHz, z]                            -> Association of fitted params
     FitCPE[freqHz, z, opts...]
   Options:
     "Weighting"     -> "Modulus" (default, 1/|Z|^2) | "Unit"
     "MaxIterations" -> 3000
*)

BeginPackage["CPEFit`"]

CPEImpedance::usage =
  "CPEImpedance[{Rs, Cdl, Rct, Q, n}, freqHz] returns the complex impedance of \
the Rs + (Cdl || (Rct + CPE)) circuit at each frequency (Hz).";

FitCPE::usage =
  "FitCPE[freqHz, z] fits the Rs + (Cdl || (Rct + CPE)) equivalent circuit to \
measured EIS data (frequency list in Hz, complex impedance list z) by weighted \
complex least-squares. Returns an Association with Rs, Cdl, Rct, Q, n, the \
reconstructed ZFit and fit-quality metrics (ChiSq, RelResid). Accepts options \
\"Weighting\" (\"Modulus\" | \"Unit\") and \"MaxIterations\".";

CPEImpedance2::usage =
  "CPEImpedance2[{Rs, Cdl, Rct1, Q1, n1, Rct2, Q2, n2}, freqHz] returns the \
complex impedance of the two-relaxation circuit Rs + (Cdl || (Rct1+CPE1) || \
(Rct2+CPE2)).";

FitCPE2::usage =
  "FitCPE2[freqHz, z] fits the two-faradaic-branch circuit Rs + (Cdl || \
(Rct1+CPE1) || (Rct2+CPE2)) for model discrimination against the single-branch \
FitCPE. Seeds from the single-branch fit, split into a slow and a fast branch. \
Returns Rs, Cdl, (Rct1,Q1,n1,ket1), (Rct2,Q2,n2,ket2) ordered slow->fast, plus \
ZFit/ChiSq/RelResid. ket = 1/(Q Rct).";

Begin["`Private`"]

(* ------------------------------------------------------------------ *)
(* Forward model                                                       *)
(* ------------------------------------------------------------------ *)

cpeZ[rs_, cdl_, rct_, q_, nn_, freqHz_List] := Module[{w},
  w = 2. Pi freqHz;
  rs + 1. / (I w cdl + 1. / (rct + 1. / (q (I w)^nn)))
]

CPEImpedance[{rs_, cdl_, rct_, q_, nn_}, freqHz_List] := cpeZ[rs, cdl, rct, q, nn, freqHz]

(* ------------------------------------------------------------------ *)
(* Weighting                                                           *)
(* ------------------------------------------------------------------ *)

weightVec["Unit", zd_]    := ConstantArray[1., Length[zd]]
weightVec["Modulus", zd_] := 1. / Abs[zd]^2
weightVec[_, zd_]         := 1. / Abs[zd]^2

(* ------------------------------------------------------------------ *)
(* Fit                                                                 *)
(* ------------------------------------------------------------------ *)

Options[FitCPE] = {"Weighting" -> "Modulus", "MaxIterations" -> 3000};

FitCPE[freqHz_List, z_List, OptionsPattern[]] := Module[
  {zd, w, wt, iHi, rs0, cdl0, rct0, q0, n0, obj, lrs, lcdl, lrct, lq, nn,
   findSol, result, best},

  zd = N[z];
  w  = 2. Pi N[freqHz];
  wt = weightVec[OptionValue["Weighting"], zd];

  (* ---- data-driven initial guesses ---- *)
  iHi  = First @ Ordering[freqHz, -1];                   (* highest-frequency point *)
  rs0  = Max[Min[Re[zd]], 1.*^-2];                       (* high-freq real intercept *)
  cdl0 = Min[Max[1. / (w[[iHi]] Max[Abs[Im[zd[[iHi]]]], 1.*^-6]), 1.*^-9], 1.*^-3];
  rct0 = Max[Max[Re[zd]] - rs0, 10.];                    (* span of Re[Z] *)
  q0   = cdl0;
  n0   = 0.85;

  (* ---- weighted complex objective, log10-parameters (n stays linear) ---- *)
  obj[a_?NumericQ, b_?NumericQ, c_?NumericQ, d_?NumericQ, e_?NumericQ] := Module[{zm},
    zm = cpeZ[10.^a, 10.^b, 10.^c, 10.^d, e, freqHz];
    Total[wt (Re[zm - zd]^2 + Im[zm - zd]^2)]
  ];

  (* ---- helper: package a solution + diagnostics, or $Failed ---- *)
  result[sol_] := If[!MatchQ[sol, {_?NumericQ, {__Rule}}], $Failed,
    Module[{p, rs, cdl, rct, q, nv, zfit},
      p   = sol[[2]];
      rs  = 10.^(lrs  /. p);  cdl = 10.^(lcdl /. p);  rct = 10.^(lrct /. p);
      q   = 10.^(lq   /. p);  nv  = nn /. p;
      zfit = cpeZ[rs, cdl, rct, q, nv, freqHz];
      <|"Rs" -> rs, "Cdl" -> cdl, "Rct" -> rct, "Q" -> q, "n" -> nv,
        "FreqHz" -> freqHz, "ZData" -> zd, "ZFit" -> zfit,
        "ChiSq" -> sol[[1]], "RelResid" -> Sqrt[Mean[Abs[(zfit - zd) / zd]^2]]|>
    ]
  ];

  (* ---- local fit from data-driven seeds ---- *)
  findSol = result @ Quiet @ Check[
    FindMinimum[
      {obj[lrs, lcdl, lrct, lq, nn], 0.3 <= nn <= 1.},
      {{lrs, Log10[rs0]}, {lcdl, Log10[cdl0]}, {lrct, Log10[rct0]},
       {lq, Log10[q0]}, {nn, n0}},
      MaxIterations -> OptionValue["MaxIterations"]],
    $Failed];

  (* ---- global fallback if the local fit failed or fitted poorly ---- *)
  best = If[findSol =!= $Failed && findSol["RelResid"] < 0.1,
    findSol,
    Module[{nmin},
      nmin = result @ Quiet @ Check[
        NMinimize[
          {obj[lrs, lcdl, lrct, lq, nn],
           Log10[1.*^-3] <= lrs  <= Log10[1.*^6]  &&
           Log10[1.*^-10] <= lcdl <= Log10[1.*^-2] &&
           Log10[1.]      <= lrct <= Log10[1.*^9]  &&
           Log10[1.*^-10] <= lq   <= Log10[1.*^-2] &&
           0.3 <= nn <= 1.},
          {lrs, lcdl, lrct, lq, nn},
          Method -> "DifferentialEvolution", MaxIterations -> 200],
        $Failed];
      (* keep whichever of the two converged better *)
      Which[
        nmin === $Failed, findSol,
        findSol === $Failed, nmin,
        nmin["RelResid"] < findSol["RelResid"], nmin,
        True, findSol]
    ]
  ];

  If[best === $Failed,
    Failure["FitCPE", <|"Message" -> "Fit did not converge"|>],
    best]
]

(* ------------------------------------------------------------------ *)
(* Two-relaxation variant for model discrimination:                    *)
(*   Rs + [ Cdl || (Rct1 + CPE1) || (Rct2 + CPE2) ]                     *)
(* Does a second faradaic branch fit dramatically better, resolving    *)
(* two distinct k_et?  Seeds from the single-branch fit, split slow/fast. *)
(* ------------------------------------------------------------------ *)

cpeZ2[rs_, cdl_, rct1_, q1_, n1_, rct2_, q2_, n2_, freqHz_List] := Module[{w},
  w = 2. Pi freqHz;
  rs + 1. / (I w cdl
     + 1. / (rct1 + 1. / (q1 (I w)^n1))
     + 1. / (rct2 + 1. / (q2 (I w)^n2)))
]

CPEImpedance2[{rs_, cdl_, rct1_, q1_, n1_, rct2_, q2_, n2_}, freqHz_List] :=
  cpeZ2[rs, cdl, rct1, q1, n1, rct2, q2, n2, freqHz]

cpeBranch[lr_, lq_, nn_] := Module[{r, q},
  r = 10.^lr; q = 10.^lq;
  <|"Rct" -> r, "Q" -> q, "n" -> nn, "ket" -> 1. / (q r)|>
]

Options[FitCPE2] = {"Weighting" -> "Modulus", "MaxIterations" -> 4000};

FitCPE2[freqHz_List, z_List, OptionsPattern[]] := Module[
  {zd, wt, single, rs0, cdl0, rct0, q0, n0, obj2,
   lrs, lcdl, lr1, lq1, n1, lr2, lq2, n2, pack, sol, best},

  zd = N[z];
  wt = weightVec[OptionValue["Weighting"], zd];

  (* seed from the single-branch fit, then split into a slow + a fast branch *)
  single = FitCPE[freqHz, z, "Weighting" -> OptionValue["Weighting"]];
  {rs0, cdl0, rct0, q0, n0} = If[FailureQ[single],
    {Max[Min[Re[zd]], 1.], 1.*^-7, 1.*^4, 1.*^-6, 0.85},
    {single["Rs"], single["Cdl"], single["Rct"], single["Q"], single["n"]}];

  obj2[a_?NumericQ, b_?NumericQ, c_?NumericQ, d_?NumericQ, e_?NumericQ,
       f_?NumericQ, g_?NumericQ, h_?NumericQ] := Module[{zm},
    zm = cpeZ2[10.^a, 10.^b, 10.^c, 10.^d, e, 10.^f, 10.^g, h, freqHz];
    Total[wt (Re[zm - zd]^2 + Im[zm - zd]^2)]
  ];

  pack[s_] := If[!MatchQ[s, {_?NumericQ, {__Rule}}], $Failed,
    Module[{p, rs, cdl, b, zfit},
      p   = s[[2]];
      rs  = 10.^(lrs /. p);  cdl = 10.^(lcdl /. p);
      (* order the two branches slow -> fast (ket ascending) *)
      b   = SortBy[{cpeBranch[lr1 /. p, lq1 /. p, n1 /. p],
                    cpeBranch[lr2 /. p, lq2 /. p, n2 /. p]}, #["ket"] &];
      zfit = cpeZ2[rs, cdl, b[[1]]["Rct"], b[[1]]["Q"], b[[1]]["n"],
                            b[[2]]["Rct"], b[[2]]["Q"], b[[2]]["n"], freqHz];
      <|"Rs" -> rs, "Cdl" -> cdl,
        "Rct1" -> b[[1]]["Rct"], "Q1" -> b[[1]]["Q"], "n1" -> b[[1]]["n"], "ket1" -> b[[1]]["ket"],
        "Rct2" -> b[[2]]["Rct"], "Q2" -> b[[2]]["Q"], "n2" -> b[[2]]["n"], "ket2" -> b[[2]]["ket"],
        "FreqHz" -> freqHz, "ZData" -> zd, "ZFit" -> zfit,
        "ChiSq" -> s[[1]], "RelResid" -> Sqrt[Mean[Abs[(zfit - zd) / zd]^2]]|>
    ]
  ];

  sol = pack @ Quiet @ Check[
    FindMinimum[
      {obj2[lrs, lcdl, lr1, lq1, n1, lr2, lq2, n2], 0.3 <= n1 <= 1. && 0.3 <= n2 <= 1.},
      {{lrs, Log10[rs0]}, {lcdl, Log10[cdl0]},
       {lr1, Log10[rct0]}, {lq1, Log10[q0 4.]}, {n1, n0},     (* slower branch seed *)
       {lr2, Log10[rct0]}, {lq2, Log10[q0 / 4.]}, {n2, n0}},  (* faster branch seed *)
      MaxIterations -> OptionValue["MaxIterations"]],
    $Failed];

  best = If[sol =!= $Failed && sol["RelResid"] < 0.1, sol,
    Module[{nm},
      nm = pack @ Quiet @ Check[
        NMinimize[
          {obj2[lrs, lcdl, lr1, lq1, n1, lr2, lq2, n2],
           Log10[1.*^-3] <= lrs  <= Log10[1.*^6]  && Log10[1.*^-10] <= lcdl <= Log10[1.*^-2] &&
           Log10[1.]     <= lr1  <= Log10[1.*^9]  && Log10[1.*^-10] <= lq1  <= Log10[1.*^-2] && 0.3 <= n1 <= 1. &&
           Log10[1.]     <= lr2  <= Log10[1.*^9]  && Log10[1.*^-10] <= lq2  <= Log10[1.*^-2] && 0.3 <= n2 <= 1.},
          {lrs, lcdl, lr1, lq1, n1, lr2, lq2, n2},
          Method -> "DifferentialEvolution", MaxIterations -> 250],
        $Failed];
      Which[nm === $Failed, sol, sol === $Failed, nm,
        nm["RelResid"] < sol["RelResid"], nm, True, sol]]
  ];

  If[best === $Failed, Failure["FitCPE2", <|"Message" -> "Fit did not converge"|>], best]
]

End[]
EndPackage[]
