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

Options[PlotDCTHeatmap] = {Theme -> Automatic};

PlotDCTHeatmap[goodSpecs_List, nK_Integer : 150, opts : OptionsPattern[]] := Module[
  {th, fg, plotOpts, nExp, kAllMin, kAllMax, kGrid, tHr, gGrid, pts, kMin, kMax},

  th       = ThemeFromOpts[{opts}];
  fg       = th["Fg"];
  plotOpts = FilterRules[{opts}, Options[ListDensityPlot]];

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
    plotOpts,
    ScalingFunctions   -> {None, "Log10"},
    InterpolationOrder -> 0,
    Sequence @@ ThemeChrome[th, 14, 1.2],
    FrameTicks         -> {{ThemeLogTicks[kMin, kMax, fg], None},
                           {ThemeLinTicks[Min[tHr], Max[tHr], fg], None}},
    FrameTicksStyle    -> Directive[fg, 13],
    FrameLabel         -> {
      Style["Time (hours)", 15, fg],
      Style[Row[{"k (", Superscript["s", -1], ")"}], 15, fg]
    },
    PlotLegends -> Placed[
      BarLegend[Automatic,
        LegendLabel      -> Style["g (F/decade)", 12, fg],
        LabelStyle       -> Directive[12, fg],
        LegendMarkerSize -> 200,
        LegendFunction   -> (Framed[#, Background -> th["Bg"],
          FrameStyle -> Directive[fg, AbsoluteThickness[0.5]]] &)
      ],
      Right
    ],
    PlotRangePadding -> Scaled[0.02],
    ImageSize        -> 750
  ]
]

End[]
EndPackage[]
