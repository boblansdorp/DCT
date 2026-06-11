(* ::Package:: *)

(*  DCTPlots`PlotTheme  —  shared plot theming.

    Single source of truth for colours/chrome across every DCTPlots figure.
    Plot functions read the active theme instead of hardcoding White/dark.

    Public API:
      $DCTTheme          global default theme: "Dark" (screen) or "Publication"
      Theme              per-call option on plot functions; Automatic -> $DCTTheme
      DCTThemeData[name] Association of theme colours for one theme
      ResolveTheme[x]    x = Automatic | name string | Association -> Association

    Each theme Association carries:
      "Bg"      background colour
      "Fg"      foreground colour (frame, ticks, axis labels)
      "Font"    label font family

    NOTE: every plot must set PlotTheme -> "Default" to escape the kernel-global
    $PlotTheme (e.g. "BlackBackground"), which would otherwise force white
    frames/ticks onto a Publication (white-background) figure. Use ThemeChrome
    below, which includes it.

    Ticks (ThemeLogTicks / ThemeLinTicks) are generated here directly — small,
    self-contained, no external paclet. Positions are VALUE coordinates, so they
    drop straight onto ScalingFunctions -> "Log10" and linear frame axes.
*)

BeginPackage["DCTPlots`"]

$DCTTheme::usage =
  "$DCTTheme is the active plot theme read by DCTPlots functions: \
\"Dark\" (default, for screen) or \"Publication\" (white background). \
Override a single call with the Theme option.";

Theme::usage =
  "Theme is an option for DCTPlots plotting functions. Values: \"Dark\", \
\"Publication\", or Automatic (use $DCTTheme).";

DCTThemeData::usage =
  "DCTThemeData[\"Dark\"|\"Publication\"] returns an Association of theme \
colours: \"Bg\", \"Fg\", \"Font\".";

ResolveTheme::usage =
  "ResolveTheme[Automatic] returns the $DCTTheme Association; \
ResolveTheme[name] looks up a named theme; ResolveTheme[assoc] passes it through.";

ThemeFromOpts::usage =
  "ThemeFromOpts[{opts}] extracts the Theme option (default Automatic) from an \
option list and returns the resolved theme Association. Avoids OptionValue so \
pass-through plot options (ImageSize, etc.) do not trigger OptionValue::nodef.";

ThemeChrome::usage =
  "ThemeChrome[theme] returns the common ListLinePlot chrome options for a \
theme Association (PlotTheme, Background, FrameStyle, LabelStyle, Frame, Axes). \
ThemeChrome[theme, fontSize, thickness] sets label size and frame thickness.";

ThemeLogTicks::usage =
  "ThemeLogTicks[fmin, fmax] / ThemeLogTicks[fmin, fmax, color] returns outward \
log-axis FrameTicks (10^n labels, value coordinates) coloured for the theme \
(default Black). For ScalingFunctions -> \"Log10\" axes.";

ThemeLinTicks::usage =
  "ThemeLinTicks[vmin, vmax] / ThemeLinTicks[vmin, vmax, color] returns outward \
linear-axis FrameTicks (value coordinates) coloured for the theme (default Black).";

Begin["`Private`"]

If[!ValueQ[$DCTTheme], $DCTTheme = "Dark"];

DCTThemeData["Dark"] = <|
  "Bg"   -> GrayLevel[0.12],
  "Fg"   -> White,
  "Font" -> "Arial"
|>;

DCTThemeData["Publication"] = <|
  "Bg"   -> White,
  "Fg"   -> Black,
  "Font" -> "Arial"
|>;

DCTThemeData[other_] := (
  Message[DCTThemeData::badtheme, other];
  DCTThemeData["Dark"]
);
DCTThemeData::badtheme =
  "Unknown theme `1`; falling back to \"Dark\". Use \"Dark\" or \"Publication\".";

ResolveTheme[Automatic]          := DCTThemeData[$DCTTheme]
ResolveTheme[name_String]        := DCTThemeData[name]
ResolveTheme[assoc_Association]  := assoc

ThemeFromOpts[opts_List] :=
  ResolveTheme[Theme /. Flatten[opts] /. {Theme -> Automatic}]

ThemeChrome[th_Association, fontSize_ : 14, thickness_ : 1.2] := {
  PlotTheme  -> "Default",
  Background  -> th["Bg"],
  FrameStyle  -> Directive[th["Fg"], AbsoluteThickness[thickness]],
  LabelStyle  -> Directive[th["Fg"], fontSize, FontFamily -> th["Font"]],
  Frame       -> True,
  Axes        -> False
}

(* ── Outward frame ticks (self-contained) ──────────────────────────────────
   A tick is {position, label, {0, length}, {color}}.  Positions are VALUE
   coordinates (right for ScalingFunctions -> "Log10" and linear axes); the
   {0, length} length spec draws the mark *outward*; the label is coloured to
   match the frame.  Major marks are longer and carry the number label, minor
   marks are shorter and blank. *)

$majorLen = 0.016;
$minorLen = 0.008;

mkTick[x_, lbl_, len_, color_] := {x, lbl, {0, len}, {color}}

(* number label for a sub-decade (1-2-5) mark on a narrow log axis *)
subDecadeLabel[x_] := Module[{e, m},
  e = Floor[Log10[x] + 10.^-9];  m = Round[x / 10.^e];
  If[-1 <= e <= 2,
    If[e >= 0, ToString[Round[x]], ToString[NumberForm[x, {2, -e}]]],
    If[m == 1, Superscript[10, e], Row[{m, "\[Times]", Superscript[10, e]}]]
  ]
]

ThemeLogTicks[fmin_?Positive, fmax_?Positive, color_ : Black] := Module[
  {lo, hi, inR, at},
  lo = Floor[Log10[fmin] + 10.^-9];
  hi = Ceiling[Log10[fmax] - 10.^-9];
  inR[x_] := fmin (1 - 10.^-6) <= x <= fmax (1 + 10.^-6);
  at[mults_] := Select[Flatten @ Table[m 10.^e, {e, lo, hi}, {m, mults}], inR];
  If[Log10[fmax] - Log10[fmin] >= 2.5,
    (* wide axis: 10^n decade labels + unlabelled 2..9 minors *)
    Join[
      mkTick[#, Style[Superscript[10, Round[Log10[#]]], color], $majorLen, color] & /@ at[{1}],
      mkTick[#, "", $minorLen, color] & /@ at[{2, 3, 4, 5, 6, 7, 8, 9}]
    ],
    (* narrow axis: 1-2-5 number labels + unlabelled 3-4-6..9 minors *)
    Join[
      mkTick[#, Style[subDecadeLabel[#], color], $majorLen, color] & /@ at[{1, 2, 5}],
      mkTick[#, "", $minorLen, color] & /@ at[{3, 4, 6, 7, 8, 9}]
    ]
  ]
]

(* nice 1-2-5 x 10^k step giving ~5 major divisions over the range *)
niceStep[raw_?Positive] := Module[{e, f},
  e = Floor[Log10[raw]];  f = raw / 10.^e;
  10.^e Which[f < 1.5, 1., f < 3.5, 2., f < 7.5, 5., True, 10.]
]

(* Compact number label: divide out a shared power of ten (cexp) when the scale
   is very small or very large, else a plain number. dec = mantissa decimals. *)
linLabel[x_, cexp_, dec_] := Module[{m},
  m = Round[N[x] / 10.^cexp, 10.^-(dec + 3)];
  m = Which[m == 0, 0, Abs[m - Round[m]] < 10.^-9, Round[m], True, NumberForm[m, {16, dec}]];
  If[cexp == 0 || m === 0, m, Row[{m, "\[Times]", Superscript[10, cexp]}]]
]

ThemeLinTicks[vmin_?NumericQ, vmax_?NumericQ, color_ : Black] := Module[
  {step, mstep, majors, minors, inR, vmag, cexp, mdec},
  If[!(vmax > vmin), Return[{}]];
  step  = niceStep[(vmax - vmin) / 5.];
  mstep = step / 5.;
  inR[x_] := vmin - step 10.^-6 <= x <= vmax + step 10.^-6;
  majors = Select[Range[Ceiling[vmin/step - 10.^-9] step,  vmax + step 10.^-6,  step],  inR];
  minors = Complement[
    Select[Range[Ceiling[vmin/mstep - 10.^-9] mstep, vmax + mstep 10.^-6, mstep], inR],
    majors];
  (* shared exponent only for extreme scales (< 1e-3 or >= 1e5) *)
  vmag = Max[Abs /@ majors];
  cexp = If[vmag > 0 && ! (10.^-3 <= vmag < 10.^5), Floor[Log10[vmag] + 10.^-9], 0];
  mdec = Max[0, -Floor[Log10[step / 10.^cexp] + 10.^-9]];
  Join[
    mkTick[#, Style[linLabel[#, cexp, mdec], color], $majorLen, color] & /@ majors,
    mkTick[#, "", $minorLen, color] & /@ minors
  ]
]

End[]
EndPackage[]
