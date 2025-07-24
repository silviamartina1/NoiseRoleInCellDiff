#This script is able to generate the bifurcation diagram in Fig. 2a of the main text 
from cmath import sqrt
import numpy as np
import sympy
import matplotlib.patches as mpatches
from pylab import *
from sympy import *
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
sympy.init_printing()
# Define variable as symbols for sympy #Python is able to symbolic calculation
X,Y = symbols('X Y ',real=True)
#parameters setting
gamma=0.005
A=0.5

#we study the one-stage Toggle Swich for two competitive molecules X and Y in the general symmetric case
K_max=A/(4*gamma)

K_max_2= 4* K_max #max value we choose for this parameter
delta=0.01 #risolution by varying K
K_values=np.arange(0,K_max,delta)
#second model ODE systems
#Here the analytical solution for X and Y concentrations related to the 5 identified steady states
def equilibrium_points_X(KX):
    X1=0.
    X2= 2*K_max*(1-KX/(2*K_max)+sqrt(1-KX/K_max))     #Y_max/2+1/2*np.sqrt(Y_max**2-4*KX*KY)
    X3= 2*K_max*(1-KX/(2*K_max)-sqrt(1-KX/K_max))     #Y_max/2-1/2*np.sqrt(Y_max**2-4*KX*KY)
    X4=0.
    X5=4*K_max-KX

    return X1,X2,X3,X4,X5

def equilibrium_points_Y(KX):
    Y1=0
    Y2= 2*K_max*(1-KX/(2*K_max)-sqrt(1-KX/K_max))  #Y_max/2-1/2*np.sqrt(Y_max**2-4*KX*KY)
    Y3= 2*K_max*(1-KX/(2*K_max)+sqrt(1-KX/K_max))  #Y_max/2+1/2*np.sqrt(Y_max**2-4*KX*KY)
    Y4=4*K_max-KX
    Y5=0

    return Y1,Y2,Y3,Y4,Y5



def nature_critical_point(KX,x0,y0):

    A = 0.5
    #symmetric case KX=KY=K
    dXdt = A*KX* X/((KX+ X)*(KX+Y)) - gamma * X

    dYdt = A*KX* Y/((KX + Y)*(KX+X)) - gamma * Y
    # Symbolic expression of the matrix
    sys = sympy.Matrix([dXdt, dYdt])  # this is the column of the differentials
    var = sympy.Matrix([X, Y])  # it the column of variables
    jac = sys.jacobian(var)  # compute the Jacobian respect to the variable X,x in a symbolic way

    jac_criticalpoint = jac.subs({X:x0,Y:y0}) #substitute at symbolic Jacobian first critical point
    A = np.array(jac_criticalpoint).astype(np.float64)
    eigenvalues= np.linalg.eigvals(A)
    # collect the eigenvalues in a array and analyze the real part of the rightmost eigenvalue
    #print(A)
    #real_part_eigenvalues = eigenvalues.real
    real_part_eigenvalues = np.round(eigenvalues.real,3)
    #print(real_part_eigenvalues)

    #now define the nature of the point according the max value of real_part_eigenvalues
    if real_part_eigenvalues[0]==0 or real_part_eigenvalues[1]==0:

        if np.max(real_part_eigenvalues)>0:
            nature=1.
        if np.max(real_part_eigenvalues)==0:
            nature=2.

    if real_part_eigenvalues[0]==0 and real_part_eigenvalues[1]==0:
        nature=3.

    elif max(real_part_eigenvalues)>0:
        nature=4.
    elif max(real_part_eigenvalues) < 0:
        nature = 5.

    return nature
    #print(real_part_eigenvalues[0],real_part_eigenvalues[2]

#np.count_nonzero(real_part_eigenvalues)<=len(real_part_eigenvalues

X0_Y=np.zeros((len(K_values),5))
X0_X=np.zeros((len(K_values),5))


for k in range(len(K_values)):
    X0_X[k]=equilibrium_points_X(K_values[k]) #collection of any K_values in the range, the corresponding value for X concentration for all 5 states
    X0_Y[k]=equilibrium_points_Y(K_values[k])

#collection of all the info inside one multi-dimensional array, remember the order of insert ,call it CP_total_nature

CP_total_nature=np.zeros((5,5,len(K_values))) #order for the states: X2,X3,X4,X5 ; order for nature: 0:stable,1:unstable,2:simplify.
CP_total_nature[:]=np.nan #initialization
for i in range(len(K_values)):
    protein_X_array=np.zeros(5)
    protein_Y_array =np.zeros(5)

    protein_X_array = equilibrium_points_X(K_values[i]) #return the 5 concentration of X for the 5 steady states
    protein_Y_array = equilibrium_points_Y(K_values[i]) #return the 5 concentration of Y for the 5 steady states
    if K_values[i]<=K_max: #(all the steady states exists )
        for j in range(5):
            Nature_K = nature_critical_point(K_values[i], protein_X_array[j], protein_Y_array[j])
            if Nature_K ==1:
                CP_total_nature[0][j][i] = protein_X_array[j]
            if Nature_K == 2:
                CP_total_nature[1][j][i] = protein_X_array[j]
            if Nature_K == 3:
                CP_total_nature[2][j][i] = protein_X_array[j]
            if Nature_K == 4:
                CP_total_nature[3][j][i] = protein_X_array[j]
            if Nature_K == 5:
                CP_total_nature[4][j][i] = protein_X_array[j]
    if K_values[i] >K_max:
        for j in(0,3,4):
            Nature_K = nature_critical_point(K_values[i],protein_X_array[j], protein_Y_array[j])
            if Nature_K == 1:
                CP_total_nature[0][j][i] = protein_X_array[j]
            if Nature_K == 2:
                CP_total_nature[1][j][i] = protein_X_array[j]
            if Nature_K == 3:
                CP_total_nature[2][j][i] = protein_X_array[j]
            if Nature_K == 4:
                CP_total_nature[3][j][i] = protein_X_array[j]
            if Nature_K == 5:
                CP_total_nature[4][j][i] = protein_X_array[j]




#for plotting Bifurcation Diagram as X or Y concentration run below

#X-concentration
mpl.rcParams.update(mpl.rcParamsDefault)
VectEquilPointsX=np.vectorize(equilibrium_points_X)
K_values=np.arange(0,K_max,delta)
#for plotting just the stable points
plt.rc('font', size=12)
plt.plot(K_values,VectEquilPointsX(K_values)[1],'hotpink',label= r"$\bar{X}_2$",linewidth=3)
plt.plot(K_values,VectEquilPointsX(K_values)[2],'orangered',label=r"$\bar{X}_3$",linewidth=3)
K_values=np.arange(K_max,K_max_2,delta)
plt.plot(K_values, VectEquilPointsX(K_values)[3],'darkgreen',label=r"$\bar{X}_4$",linewidth=3)
plt.plot(K_values,VectEquilPointsX(K_values)[4],'royalblue',label=r"$\bar{X}_5$",linewidth=3)
plt.axvline(x=K_max,ls='dashed',color='black',linewidth=2)

plt.annotate(r"$K_{max}$",xy=(25.,35),xytext=(K_max-3,-10),fontsize=15)

plt.ylabel(r'Concentration $X$',fontsize=16)
plt.title(r'$Bifurcation$ $diagram$',fontsize=15)
plt.legend()
plt.show()

#Y-concentration
mpl.rcParams.update(mpl.rcParamsDefault)
VectEquilPointsY=np.vectorize(equilibrium_points_Y)
K_values=np.arange(0,K_max+delta,delta)
#K_values=np.arange(K_max,K_max_2,delta) #here for plotting the secon stability regime for S4 and S5
#for plotting just the stable points
plt.rc('font', size=12)
plt.plot(K_values,VectEquilPointsY(K_values)[1],'pink',label= r"$S_2$",linewidth=3)
plt.plot(K_values,VectEquilPointsY(K_values)[2],'darkorange',label=r"$S_3$",linewidth=3)
K_values_2=np.arange(K_max,K_max_2,delta) #here for plotting the secon stability re
plt.plot(K_values_2, VectEquilPointsY(K_values_2)[3],'royalblue',linewidth=3)
plt.plot(K_values_2,VectEquilPointsY(K_values_2)[4],'darkgreen',linewidth=3)
plt.axvline(x=K_max,ls='dashed',color='black',linewidth=2)
plt.plot(K_values, VectEquilPointsY(K_values)[3],'royalblue',linewidth=3)
plt.plot(K_values,VectEquilPointsY(K_values)[4],'darkgreen',linewidth=3)
plt.plot(K_values_2, VectEquilPointsY(K_values_2)[3],'royalblue',linewidth=3)
plt.plot(K_values_2,VectEquilPointsY(K_values_2)[4],'darkgreen',linewidth=3)
plt.axvline(x=K_max,ls='dashed',color='black',linewidth=2)

plt.annotate(r"$K_{max}$",xy=(25.,35),xytext=(K_max-3,-15),fontsize=22)
plt.rc('font', size=22)
plt.xlabel(r'$K$',fontsize=22)
plt.ylabel(r'Concentration $Y$',fontsize=16)
plt.title(r'$Bifurcation$ $diagram$',fontsize=15)
plt.legend(fontsize=22)
plt.show()
plt.savefig('bifurcation_diagram_Fig2a.svg', format='svg')
plt.savefig('bifurcation_diagram_Fig2a.png', format='png')


#for plotting the nature of the Steady States (Additional analysis NOT present in the paper)
mpl.rcParams.update(mpl.rcParamsDefault)
K_values=np.arange(0,K_max,delta)
for i in range(5):
    plt.plot(K_values, CP_total_nature[0][i] , 'black', linewidth=3)
    plt.plot(K_values, CP_total_nature[1][i], 'grey', linewidth=3)
    plt.plot(K_values, CP_total_nature[2][i], 'red', linewidth=3)
    plt.plot(K_values, CP_total_nature[3][i], 'yellow', linewidth=3)
    plt.plot(K_values, CP_total_nature[4][i], 'orange', linewidth=3)

plt.plot(K_values, CP_total_nature[3][0], 'orange', linewidth=3) #questo e relativo al punto X0


#plt.plot(K_values, CP_total_nature[0][1] , 'black', linewidth=3)
#plt.plot(K_values, CP_total_nature[1][1], 'grey', linewidth=3)
#plt.plot(K_values, CP_total_nature[2][1], 'red', linewidth=3)

#for i in K_values:
    #nature_critical_point(i, equilibrium_points_X(i)[2],equilibrium_points_Y(i)[2])


plt.axvline(x=25,ls='dashed',color='black',linewidth=2)

plt.annotate(r'$K_{max}$',xy=(25.,35),xytext=(K_max-3,-12),fontsize=15)


plt.title(r'Steady states nature',fontsize=15)

plt.xlabel(r'$K$',fontsize=15)
plt.ylabel("X-Concentration ",fontsize=15)

black_patch = mpatches.Patch(color='black', label=r'Unstable Points (Not Asymptotically)')
grey_patch = mpatches.Patch(color='grey', label=r'Stable Points (Not Asymptotically)')
red_patch = mpatches.Patch(color='red', label=r' Null Cycle')
orange_patch = mpatches.Patch(color='orange', label=r' Unstable Points')
plt.legend(handles=[grey_patch,red_patch,orange_patch])

#plt.title("Bifurcation Diagram", pad=30,fontsize=20)

plt.show()
