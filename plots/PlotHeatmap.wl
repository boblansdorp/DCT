(* ::Package:: *)

(* DCTPlots`PlotHeatmap
   2-D heatmap of g(k) over time: x = time (hours), y = k (s^-1), colour = g.
   Interpolates each spectrum onto a common log-k grid.

   Public API:
     PlotDCTHeatmap[goodSpecs]
       -> ListDensityPlot with log k axis and time in hours on x

   Default colour scheme is dark (white-on-dark) for wolfbook display.
   For light-background export pass: Background->White,
     FrameStyle->Directive[Black,AbsoluteThickness[1.2]],
     LabelStyle->Directive[Black,14,FontFamily->"Arial"]
*)

BeginPackage["DCTPlots`"]

PlotDCTHeatmap::usage =
  "PlotDCTHeatmap[goodSpecs] returns a 2-D density plot of g(k) vs time. \
x-axis: time (hours), y-axis: k (s^-1, log scale), colour: g (F/decade)."

Begin["`Private`"]

$darkBg = GrayLevel[0.12];

PlotDCTHeatmap[goodSpecs_List, nK_ : 150, opts : OptionsPattern[]] := Module[
  {nExp, kAllMin, kAllMax, kGrid, tHr, gGrid, pts, kMin, kMax},

  nExp = Length[goodSpecs];
  If[nExp == 0, Return[$Failed]];

  kAllMin = Min[Flatten[(1. / #Spec["Tau"] &) /@ goodSpecs]];
  kAllMax = Max[Flatten[(1. / #Spec["Tau"] &) /@ goodSpecs]];

  kGrid = 10.^Subdivide[Log10[kAllMin], Log10[kAllMax], nK - 1];

  tHr = (#Spec["FinishTimeS"] / 3600. &) /@ goodSpecs;

  gGrid = Table[
    Module[{spec, kData, gData, ord, interp},
      spec  = goodSpecs[[i, "Spec"]];
      kData = 1. / spec["Tau"];
      gData = spec["g"];
      ord   = Ordering[kData];
      interp = Interpolation[
        Transpose[{Log10[kData[[ord]]], gData[[ord]]}],
        InterpolationOrder -> 1
      ];
      interp /@ Log10[kGrid]
    ],
    {i, nExp}
  ];

  pts = Flatten[
    Table[{tHr[[i]], kGrid[[j]], gGrid[[i, j]]}, {i, nExp}, {j, nK}],
    1
  ];

  {kMin, kMax} = MinMax[kGrid];

  ListDensityPlot[
    pts,
    ScalingFunctions   -> {None, "Log10"},
    InterpolationOrder -> 0,
    Frame              -> True,
    Axes               -> False,
    Background         -> $darkBg,
    FrameStyle         -> Directive[White, AbsoluteThickness[1.2]],
    LabelStyle         -> Directive[White, 14, FontFamily -> "Arial"],
    FrameLabel         -> {
      Style["Time (hours)", 15, White],
      Style[Row[{"k (", Superscript["s", -1], ")"}], 15, White]
    },
    PlotLegends -> Placed[
      BarLegend[Automatic,
        LegendLabel      -> Style["g (F/decade)", 12, White],
        LabelStyle       -> Directive[12, White],
        LegendMarkerSize -> 200,
        LegendFunction   -> (Framed[#, Background -> $darkBg,
          FrameStyle -> Directive[White, AbsoluteThickness[0.5]]] &)
      ],
      Right
    ],
    PlotRangePadding -> Scaled[0.02],
    ImageSize        -> 750,
    opts
  ]
]

End[]
EndPackage[]
