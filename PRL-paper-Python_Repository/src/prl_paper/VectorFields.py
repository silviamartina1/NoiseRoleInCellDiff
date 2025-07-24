#this script is for generating vector fields (FIGs.2a, S1a,S1b,S1c,S1d,S6a,S6b,S6c)

import numpy as np
import sympy
from pylab import *
from sympy import *
import matplotlib.pyplot as plt  # Plotting library
from sympy import *
from pylab import *
#import pylab


sympy.init_printing()
# Define variable as symbols for sympy #Python is able to to symbolic calculation
X, Y = symbols('X Y', real=True)
# parameters setting
alpha = 0.5
gamma = 0.005
Kmax = alpha / (4 * gamma)
# Kmax=alpha/(gamma)
epsilon = 0
delta_x = 0. # =0 to generate the figures in the paper.
delta_y = 0. # =0 to generate the figures in the paper.
KX = 80
KY =80
K = 50
# XSS=2*Kmax*((1-K/(2*Kmax))+np.sqrt(1-K/Kmax))
# YSS=2*Kmax*((1-K/(2*Kmax))-np.sqrt(1-K/Kmax))
# second model ODE systems Mono-stsge Toggle-Switch

# MODEL WITH BASELINE OF TRASCRIPTION AS TRANSLATION EFFECT
# dXdt =alpha * KY/(KY+Y-delta_y)*(X-delta_x)/(KX+X-delta_x) -gamma *(X-delta_x)
# dYdt= alpha * KX/(KX+X-delta_x)*(Y-delta_y)/(KY+Y-delta_y)- gamma *(Y-delta_y)

# MODEL WITH BASELINE OF TRASCRIPTION WHICH MODIFIES THE DYNAMICS
dXdt = alpha * KY / (KY + Y) * X / (KX + X) - gamma * X + delta_x
dYdt = alpha * KX / (KX + X) * Y / (KY + Y) - gamma * Y + delta_y

# FIND THE CRITICAL POINTS associated to the system

criticalpoints = sympy.nonlinsolve([dXdt, dYdt], [X, Y])

# Symbolic expression of the matrix
sys = sympy.Matrix([dXdt, dYdt])  # this is the column of the differentials
var = sympy.Matrix([X, Y])  # it the column of variables
jac = sys.jacobian(var)  # compute the jacobian respect to the variable X,x

# You can convert jac to a function:
# jacobian_ODE_crossinhibition= sympy.lambdify((X, Y, KX, KY), jac, dummify=False)


# END OF SYMBOLIC PART

# Build the grid
Y_values ,X_values =np.mgrid[-5:100:.01, -5:100:.01]
#X_values ,Y_values =np.meshgrid(np.arange(0,100,0.01),np.arange(0,100,0.01)) #this range has be chosen according the position of the critical point X2 and X3
#X_values, Y_values = np.meshgrid(np.arange(0, 60, 0.01), np.arange(0, 60, 0.01))  # this the range for X4 and X5

#Y_values ,X_values =np.mgrid[0:30:.01, 0:30:.01]
###TRAJECTORIES AROUND CRITICAL POINT
# FIRST MODEL
Xdot = alpha * KY / (KY + Y_values - delta_y) * (X_values - delta_x) / (KX + X_values - delta_x) - gamma * (
            X_values - delta_x)
Ydot = alpha * KX / (KX + X_values - delta_x) * (Y_values - delta_y) / (KY + Y_values - delta_y) - gamma * (
            Y_values - delta_y)
# SECOND MODEL
# Xdot=alpha * KY/(KY+Y_values)*X_values/(KX+X_values) -gamma *(X_values)+delta_x
# Ydot=alpha * KX/(KX+X_values)*Y_values/(KY+Y_values)- gamma *(Y_values)+delta_y


speed = np.sqrt(Xdot ** 2 + Ydot ** 2)
fig, ax = plt.subplots()
plt.rcParams["figure.figsize"] = (8., 8.)
lw = 3 * speed / speed.max()
# mtp.pyplot.streamplot(X_values,Y_values,Xdot,Ydot,density=[3,3],color='grey')#,linewidth=lw)
plt.streamplot(X_values, Y_values, Xdot, Ydot, density=3, color='grey')
# plt.plot(criticalpoints.args[2][0],criticalpoints.args[2][1],'ro')
# here for the two branches of the stability manifold K=25
Kplot=np.arange(0,25.01,0.01)
#Kplot=np.arange(25,K,0.01)
# ZerosPlot=np.zeros(len(Kplot))
plt.plot(2*Kmax*((1-Kplot/(2*Kmax))+np.sqrt(1-Kplot/Kmax)),2*Kmax*((1-Kplot/(2*Kmax))-np.sqrt(1-Kplot/Kmax)),'o',markersize=3,color='pink')
plt.plot(2*Kmax*((1-Kplot/(2*Kmax))-np.sqrt(1-Kplot/Kmax)),2*Kmax*((1-Kplot/(2*Kmax))+np.sqrt(1-Kplot/Kmax)),'o',markersize=3,color='orange')
# plt.plot(4*Kmax-Kplot,4*Kmax-Kplot)
# K=Kmax;
XSS=2*Kmax*((1-K/(2*Kmax))+np.sqrt(1-K/Kmax))
YSS=2*Kmax*((1-K/(2*Kmax))-np.sqrt(1-K/Kmax))

# here the dashed line for the stability manifold vector field
# Kplot2=np.arange(0,25,0.01)
# plt.plot(2*Kmax*((1-Kplot2/(2*Kmax))+np.sqrt(1-Kplot2/Kmax)),2*Kmax*((1-Kplot2/(2*Kmax))-np.sqrt(1-Kplot2/Kmax)),'--',markersize=0.1,color='k')
# plt.plot(2*Kmax*((1-Kplot2/(2*Kmax))-np.sqrt(1-Kplot2/Kmax)),2*Kmax*((1-Kplot2/(2*Kmax))+np.sqrt(1-Kplot2/Kmax)),'--',markersize=0.1,color='k')

# it is not a linear relation!
# Xplot=np.arange(0,4*Kmax-K,0.01)
# plt.plot(Xplot,4*Kmax-K-Xplot,'.',markersize=1.5,color='k')
plt.plot(XSS,YSS,'o',markersize=10,color="hotpink")
plt.plot(YSS,XSS,'o',markersize=10,color="orangered")
plt.annotate(r'$S_2$', (XSS,YSS+4), fontsize=25)
plt.annotate(r'$S_3$', (YSS+4,XSS),fontsize=25, weight='bold')
plt.plot(0, 4 * Kmax - KY, 'o', markersize=10, color="blue")
plt.plot(4 * Kmax - KX, 0, 'o', markersize=10, color="green")
plt.plot(0, 0, 'o', markersize=10, color="black")
plt.annotate(r'$S_5$',( 4 * Kmax - KX+2, 2), fontsize=25)
plt.annotate(r'$S_4$', (1, 4 * Kmax - KY+2), fontsize=25)
plt.annotate(r'$S_1$', (1, 1), fontsize=23, weight='bold')
# here when K=Kmax
# plt.annotate(r'$S_2 \equiv S_3 $', (XSS,YSS+2),fontsize=25) #K=K_max
# plt.plot(YSS,XSS,'o',markersize=10,color="orangered")
# X_line=np.arange(0,100,0.01)
# plt.plot(X_line ,4*Kmax-K-X_line,'--',color='black')
# ylim(bottom=-1)

plt.xlabel(r'$X$ Concentration', fontsize=20)
plt.ylabel(r'$Y$ Concentration', fontsize=20)
#plt.title(r'$K=50$', fontsize=25)
plt.title(r'$K_{X}=80, K_{Y}=80$', fontsize=20)
plt.xticks(np.arange(0, 100, step=20),fontsize=22)
plt.yticks(np.arange(0, 100, step=20),fontsize=22)
plt.rcParams["figure.figsize"] = (8., 8.)
plt.show()
plt.savefig('Vector_fields.jpeg',dpi=400)
# here for anaylsising dynamics with baseline
# stable_state=np.roots([gamma,2*gamma*K-delta_x,gamma*K**2-2*K*delta_x-alpha*K,-delta_x*K**2])
# plt.plot(stable_state[1],stable_state[1],'o',color="blue",markersize=10)

# stable_state=np.roots()
#here below to explore the effect of adding the baseline
plt.plot(Kmax - KX, delta_y, 'o', color="darkgreen", markersize=10)
plt.plot(delta_x, Kmax - KY, 'o', color="royalblue", markersize=10)
plt.annotate(r'$S_5$', (4 * Kmax - KX + 2, -1), fontsize=20, weight='bold')
plt.annotate(r'$S_4$', (2, 4 * Kmax - KY - 1), fontsize=20)
plt.xlabel(r'$X$ Concentration', fontsize=20)
plt.ylabel(r'$Y$ Concentration', fontsize=20)
plt.title(r'$K_{X}=50, K_{Y}=80,  \delta_x=\delta_y=0.01$', fontsize=20)
matplotlib.rc('xtick', labelsize=20)
matplotlib.rc('ytick', labelsize=20)
plt.rcParams["figure.figsize"] = (8, 5.5)
plt.show()
fig.savefig('KX_50_KY_80_delta_x_0.01.pdf')
