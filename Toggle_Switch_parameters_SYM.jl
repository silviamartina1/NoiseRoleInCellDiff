#parameters setting for SYMmetric case (on the stable manifold)
A=Float64(0.5); #transcription rate
gamma=Float64(0.005); #protein degradation 5 * 10^-3 to 5 *10^-6
Kmax=A/(4*gamma);
K=25;
KX=K;
KY=K;
#Steady state concentration depending on K
XSS=2*Kmax*(1-K/(2*Kmax)+sqrt(1-K/Kmax));
YSS=2*Kmax*(1-K/(2*Kmax)-sqrt(1-K/Kmax));
#D=0.001
D=0.0001 #for K=20
#diffusive coefficient
mu_x=mu_y=sqrt(2*D)
M=Int64(1000000); #number of cells
#Initizialization of the variables

#Matrix inizialization
Mono_TS_matrix=zeros((M,2));
X0= XSS;#25;#52.3;#36;#20.73; #20
Y0= YSS;#15;#7.6;#16;#29.67; #50
X_in=X0*ones(M);
Y_in=Y0*ones(M);
Mono_TS_matrix[:,1]=X_in;
Mono_TS_matrix[:,2]=Y_in;
