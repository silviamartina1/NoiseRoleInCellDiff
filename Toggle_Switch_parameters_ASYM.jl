#asymmetric study of variance
A=Float64(0.5); #protein degradation 5 * 10^-3 to 5 *10^-6

#MRNA LEVEL

gamma=Float64(0.005); #transcription rate

Kmax=Float64(A/(gamma));
KX=30;
KY=30;

delta_x=0.01;
delta_y=0.01;
YSS=A*KX/(gamma*(KX+delta_x))-KY+delta_y/gamma; #S4 is a point on the stability manifold by varying K
XSS=delta_x;
#Matrix inizialization
#X0=XSS;
#Y0=YSS;
#X0=15.7;#initial point with the baseline
#Y0=15.7;
M=Int64(1000000); #number of trajectories

D=0.001
alpha=acos(1/(sqrt(1+((Kmax-KX)/(Kmax-KY))^2)));
#diffusive coefficient
mu_x=mu_y=sqrt(2*D)


#X_in=X0*ones(M);
#Y_in=Y0*ones(M);
#Matrix inizialization
#Mono_TS_matrix=zeros((M,2));

#Mono_TS_matrix[:,1]=X_in;
#Mono_TS_matrix[:,2]=Y_in;

#here the parameters to generate FIG. S7
KX=80
KY=80
K=80
X0=0
Y0=Kmax-K
X_in=X0*ones(M);
Y_in=Y0*ones(M);
#Matrix inizialization
Mono_TS_matrix=zeros((M,2));

Mono_TS_matrix[:,1]=X_in;
Mono_TS_matrix[:,2]=Y_in;
