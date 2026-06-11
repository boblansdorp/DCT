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
log-axis FrameTicks (CustomTicks, 10^n labels) coloured for the theme (default \
Black). For ScalingFunctions -> \"Log10\" axes; positions are real coordinates.";

ThemeLinTicks::usage =
  "ThemeLinTicks[vmin, vmax] / ThemeLinTicks[vmin, vmax, color] returns outward \
linear-axis FrameTicks (CustomTicks) coloured for the theme (default Black).";

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

(* ── Outward ticks (CustomTicks backend) ───────────────────────────── *)
(* CustomTicks colours the tick *marks* via MajorTickStyle/MinorTickStyle;
   the label text is coloured here by wrapping each non-empty label. *)

colorTicks[ticks_, color_] :=
  Replace[ticks,
    {x_, lbl_, rest___} :> {x, If[lbl === "", "", Style[lbl, color]], rest},
    {1}]

(* Decade-only ticks are too sparse on narrow (< 2.5 decade) log axes — a Bode
   panel can end up with one or even zero number labels. For such ranges we label
   the 1-2-5 sub-decade positions (3/4/6-9 as unlabelled minor marks) so the axis
   carries enough ticks and numbers. Positions are real coords (ScalingFunctions). *)
subDecadeLabel[x_] := Module[{e, m},
  e = Floor[Log10[x] + 10.^-9];  m = Round[x / 10.^e];
  If[-1 <= e <= 2,
    If[e >= 0, ToString[Round[x]], ToString[NumberForm[x, {2, -e}]]],
    If[m == 1, Superscript[10, e], Row[{m, "\[Times]", Superscript[10, e]}]]
  ]
]

subDecadeLogTicks[fmin_, fmax_, color_] := Module[{lo, hi, pick},
  lo = Floor[Log10[fmin] + 10.^-9];
  hi = Ceiling[Log10[fmax] - 10.^-9];
  pick[mults_] := Select[
    Flatten @ Table[m 10.^e, {e, lo, hi}, {m, mults}],
    fmin (1 - 10.^-6) <= # <= fmax (1 + 10.^-6) &];
  Join[
    ({#, Style[subDecadeLabel[#], color], {0, 0.016}, {color}} &) /@ pick[{1, 2, 5}],
    ({#, "", {0, 0.008}, {color}} &) /@ pick[{3, 4, 6, 7, 8, 9}]
  ]
]

ThemeLogTicks[fmin_?Positive, fmax_?Positive, color_ : Black] :=
  If[Log10[fmax] - Log10[fmin] >= 2.5,
    colorTicks[
      CustomTicks`LogTicks[fmin, fmax,
        CustomTicks`LogPlot       -> True,    (* emit real coords for ScalingFunctions *)
        CustomTicks`TickDirection -> Out,
        CustomTicks`MajorTickLength -> 0.016,
        CustomTicks`MinorTickLength -> 0.008,
        CustomTicks`MajorTickStyle  -> {color},
        CustomTicks`MinorTickStyle  -> {color}],
      color],
    subDecadeLogTicks[fmin, fmax, color]
  ]

ThemeLinTicks[vmin_?NumericQ, vmax_?NumericQ, color_ : Black] :=
  colorTicks[
    CustomTicks`LinTicks[vmin, vmax,
      CustomTicks`TickDirection -> Out,
      CustomTicks`MajorTickLength -> 0.016,
      CustomTicks`MinorTickLength -> 0.008,
      CustomTicks`MajorTickStyle  -> {color},
      CustomTicks`MinorTickStyle  -> {color}],
    color]

End[]
EndPackage[]
