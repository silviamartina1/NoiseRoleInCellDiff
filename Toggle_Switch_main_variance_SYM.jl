println("This is a  julia code to study variance on the stable manifold \n")
#this script is for analysing variance trends on the stable manifold
#to generate FIG. S4

#needed packages
using Distributed
using Plots
using DelimitedFiles
using LsqFit
using LaTeXStrings
using Base.Threads
using Distributions
using StatsBase, StatsPlots
using TickTock
using Random; Random.seed!(convert(Int,ceil(abs(time()- ceil(time()))*1000))); #or well defined Random.seed!(0);
using Plots;gr()
using DelimitedFiles
using DataFrames
using CSV
#using PlotlyJS
#using Gadfly
using BenchmarkTools

#nWorkers = 5 # set the number of workers here
#addprocs(nWorkers) #
print("Number of cores:",nprocs());
print("Number of workers:",nworkers());
@everywhere include("./Toggle_Switch_parameters_SYM.jl")
@everywhere include("./Toggle_Switch_functions.jl")

### toggle variables and parameters


### Simulation timeXo and loop variables
# has to be defined before fhn parameters for noise discretization
#Here Local variable
@everywhere dt= Float64(0.01); #time step of simulation equal to delta  t of Brownian path for variance
@everywhere T = Float64(8000); # total sim time for variance
@everywhere preT = Float64(0.); # equilibrate sim time
@everywhere t = Float64(0.0);  # current time
@everywhere ds=Float64(25); #sample rate for writing to file
#@everywhere ds=Float64(10); #sample rate for writing to file when you need just the final time
 #number of trajectories

# parameters for ploting
 #here local Variables
ni=Float64(1);  #counting variable for output
ns= Float64(1);# counting variable for output
nh=Float64(1); #counting variable for histogram output
dh= T/3;

#filename to open for printing the output
@everywhere fname = "Mono_Toggle_on_manifold_funcmodel_KX_$KX"*"_$M"*"_T_$T"*"_new.txt"
io =open(fname, "w");

@time begin
### time loop STOCHASTIC SOLUTION FOR ENSEMBLE SIMULATION
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

      #Mono_TS_matrix += mono_togswich_stoch_array(Mono_TS_matrix,A,KX,KY,gamma,dt);# calls the Toggle-Switch-ODE deterministic
      Mono_TS_matrix += mono_togswich_stoch_array(Mono_TS_matrix,A,KX,KY,gamma,dt);
      #Mono_TS_matrix= abs.(Mono_TS_matrix)




      global Mono_TS_matrix



      #ns=outputMonoTS_matrix!(ns,io, t, preT,ds) #here to collect the values for X and Y to every trajectory
      #ns=outputMonoTS_XY_bar_matrix!(ns,io, t, preT,ds)#here to collect the values for X^ and Y^ to every trajectory
      #ns=outputVariance_matrix!(ns,io, t, preT,ds)
      #ns=outputMean_matrix!(ns,io, t, preT,ds)
      #ns=outputVariance_Mean_Norm_matrix!(ns,io, t, preT,ds)
      #ns=outputVariance_Mean_matrix!(ns,io, t, preT,ds) #here for the variance along X and Y respect to X0 and Y0
      ns=outputVariance_Mean_ontheSM_matrix!(ns,io, t, preT,ds,theta)#here the variance in the new reference system tangent to stability manifold







end  #end t while loop
close(io)
end


#after filling the fname file the following code is for analysing the variance
gname ="Mono_Toggle_on_manifold_funcmodel_KX_25_1000000_T_8000.0_new.txt";
sqas=readdlm(gname);
Plots.PlotlyBackend()
#In gname the info are in columns t X Y  and then X Y  for any trajectory when
#ns=outputMonoTS_matrix!(ns,io, t, preT,ds)
size(sqas)[2]; #total number of columns in gname

#here to evaluate the variance trend when ns=outputVariance_matrix! and
#we want to plot
K=25; #change anytime gname changes

#X0=2*Kmax*(1-K/(2*Kmax)+sqrt(1-K/Kmax));
Y0=2*Kmax*(1-K/(2*Kmax)-sqrt(1-K/Kmax));
Plots.plot(sqas[:,1],sqas[:,3].-Y0,xticks = 0:15:100, ylims=(-0.015,0.015),markershape=:circle,markersize=:4,linestyle=:dot,lw=:2,linecolor =:black,markercolor=:black,ylabel=L"\Delta({\mu}_Y(t);Y_0)",xlabel=L"\textit{t}",label=L"K=5",xguidefontsize=18,
xtickfont=font(10),ytickfont=font(10),yguidefontsize=18,legend=:false)

# One-dimensional & Bidimensional histogram P(X,Y)
X_T=sqas[1,2:2:size(sqas)[2]];
Y_T=sqas[1,3:2:size(sqas)[2]]; #all the X concentration at time T defined by sqas[*,] then *=size(sqas)[1] is the final time

eta_T=(X_T-Y_T)/2;
K_var=Kmax.-eta_T.^2 ./(4*Kmax);


#Plot <x^2> from the file written directly the variance
#when gname is the file written with ns=outputVariance_Mean_matrix!(ns,io, t, preT,ds)
plot(sqas[8:size(sqas)[1],1],log.(2,sqas[8:size(sqas)[1],3]),xticks = 0:15:100, markershape=:circle,markersize=:4,markercolor=:black,linestyle=:dot,lw=:2,linecolor =:black,xlabel=L"\textit{t}",xguidefontsize=18,
xtickfont=font(10),ytickfont=font(10),yguidefontsize=18,label=L"K=5")
#general behavior sigma x and y for K
#y
#plot!(sqas[1000:end,1],D*sqas[1000:end,1].* (2 .+sqas[1000:end,1].*gamma .*(-1 .+ sqrt(1- 4*KX*gamma/A)) ),lw=1.4,linestyle=:dashdotdot,linecolor = :green, label=L"K=10 \, \textit{analytical result}",legend=(.15,.88)) #y
#x
KX=25;
#plot of theoretical predictions Eq. S4.2
plot!(sqas[1:end,1],log.(2,D*sqas[1:end,1].* (2 .-sqas[1:end,1].*gamma .*(1 .+ sqrt(1- 4*KX*gamma/A)))),linestyle=:solid,linecolor =:dodgerblue3,lw=:1.8,legendfont=font(9),
label=false) #x direction

plot!(sqas[8:end,1],log.(2,D*sqas[8:end,1].* (2 .+sqas[8:end,1].*gamma .*(-1 .+ sqrt(1- 4*KX*gamma/A)))),lw=1.8,linecolor=:black,label=false
,legend=false,legendfont=font(9)) #y direction

plot!(legend=(.80,.30),legendfont=font(12),xguidefontsize=18,
xtickfont=font(10),ytickfont=font(10),yguidefontsize=18) #y direction
xlabel!(L"t")
ylabel!(L"\log_2(\sigma^2_Y(t))") #FIGs. S4a and S4b


#in the new coordinate space X^ and Y^ on the stable manifold using ns=outputVariance_Mean_ontheSM_matrix!(ns,io, t, preT,ds,theta)
#for log normal trend
Plots.plot!(sqas[1:5:size(sqas)[1],1],log.(2,sqas[1:5:size(sqas)[1],3]),xticks = 0:1000:8000, markershape=:circle,markersize=:3,markercolor=:black,linestyle=:dot,lw=:2,linecolor =:grey,xguidefontsize=18,
xtickfont=font(10),ytickfont=font(10),yguidefontsize=18,label=L"K=25")

#for variance trend (from simulation) using ns=
Plots.plot(sqas[1:size(sqas)[1],1],sqas[1:size(sqas)[1],2],xticks = 0:500:8000, markershape=:circle,markersize=:4,markercolor=:seagreen,linestyle=:dot,lw=:2,linecolor =:grey,xguidefontsize=18,
xtickfont=font(10),ytickfont=font(10),yguidefontsize=18,label=L"K=5 ,10 ,15 ,25\quad\textrm{ from simulations}")

D=0.0001
K=20; #here the parameter describing the steady states on the stability manifold
X0=2*Kmax*(1-K/(2*Kmax)+sqrt(1-K/Kmax)); #X2 is a point on the stability manifold by varying K
Y0=2*Kmax*(1-K/(2*Kmax)-sqrt(1-K/Kmax));
theta=pi+atan(1-sqrt(Kmax/X0)-(1-sqrt(Kmax/X0))/(sqrt(1+X0/Kmax-2*sqrt(X0/Kmax)))); # it is valid for any X0 except X0=25 where the function has the discontinuity
#theta=3*pi /4; #use just for K=25

final_time_point=size(sqas)[1]
final_time_point=20
#theoretical predictions in Eq.S4.4
# sigma on X^
Plots.plot(sqas[1:final_time_point,1],D*sqas[1:final_time_point,1] .*(2 .- sqas[1:final_time_point,1].*gamma .+ sqas[1:final_time_point,1].*gamma
 .* sqrt(1-4*K*gamma/A).*cos(2*theta).+
sqas[1:final_time_point,1].*(A/K-3*gamma).*sin(2*theta)),linestyle=:solid,linecolor =:black,lw=:1.8,legendfont=font(9),xlabel=L"\textit{t}",
xtickfont=font(10),xguidefontsize=18,ytickfont=font(10),yguidefontsize=18,legend=(.75,.35),label=L"K=5 ,10 ,15 ,25\quad\textrm{ from analytics}")

# sigma on Y^
Plots.plot(sqas[1:final_time_point,1],log.(2,D*sqas[1:final_time_point,1] .*(2 .- sqas[1:final_time_point,1].*gamma .- sqas[1:final_time_point,1].*gamma
.* sqrt(1-4*K*gamma/A).*cos(2*theta).-
sqas[1:final_time_point,1].*(A/K-3*gamma).*sin(2*theta))),linestyle=:solid,linecolor =:black,lw=:1.8,xguidefontsize=18,ytickfont=font(10),xlabel=L"\textit{t}",yguidefontsize=18,legendfont=font(9),label=false)

Plots.plot!(legend=(.15,.90),legendfont=font(12))
#Plots.plot!(legend=false)
xlabel!(L"t")
ylabel!(L"\log_2(\sigma^2_{\hat{Y}}(t))") #FIGs. S4c and S4d

savefig("variance_long_time_Y_bar_M_$M"*"_T_$T"*".svg")
savefig("variance_Y_bar_K_"*"_$M"*"_$T"*".svg")
savefig("variance_X_K_"*"_$M"*"_$T"*".svg")
savefig("variance_Y_K_"*"_$M"*"_$T"*".svg")

#Plot norm of the variance from the file written directly the variance
#with ns=outputVariance_Mean_Norm_matrix!

plot(sqas[1:size(sqas)[1],1],sqas[1:size(sqas)[1],2],lw=1.8,label="variance simulation",legend=(.65,.38))
#theoretical prediction of the variance norm Eq.S4.3
plot!(sqas[1:size(sqas)[1],1],2*D*sqas[1:size(sqas)[1],1].* (2 .- gamma*sqas[1:size(sqas)[1] ,1]),lw=1.8,linecolor = :red,label="norm variance",legend=(.65,.38))
