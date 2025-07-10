println("This is a  julia code to run stochastic simulation on the stability manifold for the Toggle-Switch stochastic dynamics \n")
#this script is for generating FIG.2c on the main text and FIGs. S2 and S3.
#needed packages
using Distributed
using Base.Threads
using StatsBase, StatsPlots, Distributions
using TickTock
using Random; Random.seed!(convert(Int,ceil(abs(time()- ceil(time()))*1000))); #or well defined Random.seed!(0);
using Plots;#gr()
using DelimitedFiles
using DataFrames
using CSV
using LsqFit
using LaTeXStrings
#using Gadfly
using Statistics
using BenchmarkTools

#nWorkers = 5 # set the number of workers here
#addprocs(nWorkers) #
print("Number of cores:",nprocs());
print("Number of workers:",nworkers());
@everywhere include("./Toggle_Switch_parameters_SYM.jl")
@everywhere include("./Toggle_Switch_functions.jl")

### toggle variables and parameters


### Simulation time and loop variables
# has to be defined before fhn parameters for noise discretization
#Here Local variable
@everywhere dt= Float64(0.1);#time step of simulation equal to delta  t of Brownian path
@everywhere T = Float64(8000); # total sim time
#@everywhere T = Float64(15); #time for the approximative local solution in the variance
preT = Float64(0.); # equilibrate sim time
@everywhere t = Float64(0.0);  # current time
#ds is the rate for writing to file
ds=Float64(4000);#for histogram P(X) at final time, i don't need to save info with high rate.
# parameters for ploting
 #here local Variables
ni=Float64(1);  #counting variable for output
ns= Float64(1);# counting variable for output
nh=Float64(1); #counting variable for histogram output
dh= T/3;

#filenames
@everywhere fname = "./Stationary_probability_K_$K"*"_M_$M"*"_T_$T"*"_D_$D"*"_X0_$X0"*"_Y0_$Y0"*"_A_$A"*"_gamma_$gamma"*".txt"
#gname = "Toggle_stoch_KX_$KX"*"_KY_$KY"*"_X0_$X_rand0"*"_Y0_$Y_rand0"*"_x0_$x_rand0"*"_y0_$y_rand0"*"_M_$M"*".txt"
io =open(fname, "w");
#is= open(gname,"w");

@time begin
### time loop for STOCHASTIC simulation
while t <= T + preT
       global t = t + dt
       @everywhere t
      #hist_data_X=Float64[];
      #hist_data_Y=Float64[];
      #hist_data_x=Float64[];
      #hist_data_y=Float64[];
       global ns
       @everywhere Mono_TS_matrix
       # are empty for each temporal step
      #global nh


      #dmci_stoc_matrix=@spawnat i (mci_stoc_coup_workers(i ,dmci_stoc_matrix)) for i in workers()
      #Mono_TS_matrix += mono_togswich_array(Mono_TS_matrix,A,KX,KY,gamma,dt);#calls the deterministic Toggle-Switch-ODE
      Mono_TS_matrix += mono_togswich_stoch_array(Mono_TS_matrix,A,KX,KY,gamma,dt); #calls the stochastic Toggle-Switch-ODE



      global Mono_TS_matrix



      ns=outputMonoTS_matrix!(ns,io, t, preT,ds)







end  #end t while loop
close(io)
end
#here the code after the file is written
fname_stat="./Stationary_probability_K_20_M_1000000_T_8000.0_D_0.0001_X0_7.6393202250021055_Y0_52.3606797749979_A_0.5_gamma_0.005.txt"
fname_stat="./Stationary_probability_K_25_M_1000000_T_30.0_D_0.001_X0_25_Y0_25_A_0.5_gamma_0.005.txt"
sqas=readdlm(fname_stat);
Plots.PlotlyBackend()
# One-dimensional & Bidimensional histogram P(X,Y)
X_T=sqas[size(sqas)[1],2:2:size(sqas)[2]]; #all the X concentration at yhe final time T defined by sqas[*,] then *=size(sqas)[1] is the final time
Y_T=sqas[size(sqas)[1],3:2:size(sqas)[2]];
eta_T=(X_T-Y_T)/2; #new variable
K_var=Kmax.-eta_T.^2 ./(4*Kmax); #new variable

Der_X_T=-1 .+ (1 .- K_var./Kmax).**(-1/2);
Der_X_T2=-1 .- (1 .- K_var./Kmax).**(-1/2);
Der_X_Tot=vcat(Der_X_T,Der_X_T2);
#X_T=sqas[1,2:2:size(sqas)[2]]; #all the X concentration at time T=200
#Y_T=sqas[1,3:2:size(sqas)[2]];
Plots.histogram2d(X_T,Y_T,nbins=140 ,normed=true,ylabel=L"Y" ,xlabel=L"X",thickness_scaling = 1.5)#,xlim=(40,70),ylim=(10,45))#,aspect_ratio=:equal)
Plots.histogram2d(K_var,X_T,nbins=500,normed=true,ylabel=L"x",xlabel= L"\kappa",thickness_scaling = 1.5)
#the above line to generate the insert in Fig. 2c in the main text for Time fixed
histogram!(Der_X_Tot,nbins=50, normed=true,color=:grey,ylabel=L"P(\Delta x/\Delta \kappa)",xlabel= L"\Delta x/\Delta \kappa",label=L"T=30",thickness_scaling = 1.5)
Plots.savefig("./PRL_insert_in_figure_2c.svg")
#,xlim=(22,25.000),ylim=(10,45))#,xlim=(22,25.000),ylim=(10,45))
#here to plot the analytic curves of the stability manifold
K_values=collect(22:0.00002:25.0000)
Plots.plot!(K_values,2*Kmax.*(1 .- K_values./(2*Kmax).+sqrt.(1 .- K_values ./Kmax)),label="",linewidth=2,linecolor=:pink)
Plots.plot!(K_values,2*Kmax .*(1 .-K_values/(2*Kmax).-sqrt.(1 .-K_values./Kmax)),label="",linewidth=2,linecolor=:orange)
Plots.savefig("./PRL_figure_S3d_T_8000.svg")
