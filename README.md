### Lunar Cratering Asymmetry Estimating Code

this code is based on formulae from doi:10.1088/1674-4527/21/1/197(Lunar cratering asymmetries with high lunar orbital obliquity and inclination of the Moon).
The original version code of doi:10.1088/1674-4527/21/1/197 has been lost. This version is rewritten based on the integration_base. The example result is "lon_lat_am60.csv".
This can be plot through mathematica:

```Mathematica
Clear["Global`*"];
lldata = Import[NotebookDirectory[] <> "lon_lat_am60.csv", "csv"];
{lon, lat, flux, nvec} = Transpose[lldata[[All, {1, 2, 3, 4}]]];
nvec = nvec/flux * 19.0;
flux = flux/Max[flux];

(*plot the flux*)
data = Transpose[{Round[lat*180.0/\[Pi], 0.01],Round[lon*180.0/\[Pi], 0.01], nvec}];
geoData = GeoPosition[data];
plot = GeoContourPlot[geoData, GeoProjection -> "Mollweide", 
   Contours -> 20,
   PlotRange -> All, ColorFunctionScaling -> True, 
   ContourShading -> True,
   ColorFunction -> "TemperatureMap", GeoModel -> "Moon" , 
   OpacityFunction -> None, BoundaryStyle -> Black,
   GeoBackground -> None, 
   PlotLegends -> 
    BarLegend[{"TemperatureMap", {Min[nvec], Max[nvec]}}, 15, 
     LegendLabel -> "normal_vel", 
     LabelingFunction -> (Round[#1, 0.01] &)]];
Export[NotebookDirectory[] <> "plot_nvec.eps", plot]

(*plot the * cratering freq*)
frecra = nvec^0.987*flux;
frecra = frecra/Max[frecra];
geoData = GeoPosition @ Transpose[{Round[lat*180.0/\[Pi], 0.01], 
     Round[lon*180.0/\[Pi], 0.01], frecra}];
plot = GeoContourPlot[geoData, GeoProjection -> "Mollweide", 
   Contours -> 20,
   PlotRange -> All, ColorFunctionScaling -> True, 
   ContourShading -> True,
   ColorFunction -> "TemperatureMap", GeoModel -> "Moon" , 
   OpacityFunction -> None, BoundaryStyle -> Black,
   GeoBackground -> None, 
   PlotLegends -> 
    BarLegend[{"TemperatureMap", {Min[frecra], Max[frecra]}}, 15, 
     LegendLabel -> "flux", 
     LabelingFunction -> (Round[#1, 0.01] &)]];
Export[NotebookDirectory[] <> "plot_freq.svg", plot]
```