## Call SW_Disp
clearconsole()
cd("/Users/omenbonum/Documents/JULIA/SWF")
## Load Modelues
using DelimitedFiles, GMT
## Include functions
include("func/SW_Disp_light.jl") # Main Funtion for SW Forward Modelling

## Read  Data
Model=readdlm("data/jAK135.txt",Float64); #This is the model that we will compute dispersion for


## INPUTS
T0=10::Int64; #Shortes period (s)
TF=300::Int64; #Longest perdiod (s)
Tstep=1::Int64 # Period Step
Nmodes=3::Int64; #Number of modes (Currently 1 to 3)
## Forward modelling
@time T, SWV =SW_Disp_light(Model,T0,Tstep,TF,Nmodes);
plot(T,SWV[:,1,1],lc=:blue,R="5/405/2.5/6.2")
plot!(T,SWV[:,2,1],lc=:red)
plot!(T,SWV[:,3,1],lc=:green)
plot!(T,SWV[:,4,1],lc=:orange,show=true)
