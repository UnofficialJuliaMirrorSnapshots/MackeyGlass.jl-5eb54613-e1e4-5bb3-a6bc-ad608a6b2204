# MackeyGlass.jl
Julia scripts to generate Mackey Glass time series.

## Info
This script generates a Mackey-Glass time series using the 4th 
order Runge-Kutta method. The code is a straighforward translation 
in Julia of Matlab code, available [here](https://ww2.mathworks.cn/matlabcentral/fileexchange/24390-mackey-glass-time-series-generator?s_tid=prof_contriblnk).                                           
 

## Exemple

```julia
 using MackeyGlass
 using Plots

 T,X = MGGenerator()
 plot(T, X, label = "Mackey Glass")
```
<p align="center">
<img width="400px" src="https://github.com/JonathanCourtois/Mackey-Glass-Generator/blob/master/MGplot.png"/>
</p>
