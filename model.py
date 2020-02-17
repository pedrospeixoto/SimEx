#Micro extraction model 
# Pedro Peixoto (ppeixoto@usp.br)

import numpy as np 
import scipy as sp 
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys

import params

header='''
--------------------------------------------------------------
Micro Extraction Diffusion Model
Pedro S. Peixoto - University of São Paulo (ppeixoto@usp.br)
--------------------------------------------------------------
'''
print(header)
#----------------------------------
# Extraction mechanism parameters defined via param.py
#----------------------------------
D=params.D #diffusion coefficients
K=params.K
C=params.C
xspace=params.x
nspaces=len(D) #number of spaces
nparts=len(K)

print("You defined a mechanism with "+str(nspaces)+" space(s).")
print("Mechanism layout (x): ",xspace)
print("Initial concentrations:", C)
print("Diffusion coefficients:", D)
print("Interface coefficients:", K)

#Check dimensions
if nparts!=nspaces-1 and nspaces>1:
    print("Number of partitions must match the number of interfaces between spaces")
    print("Please re-configure parameters in params.py")
    sys.exit(-1)

#Discretize space
n = 10
dxspace = np.array([j-i for i, j in zip(xspace[:-1], xspace[1:])])
mindxspace=min(dxspace)
dx=mindxspace/n
a = xspace[0]
b = xspace[-1]
N = int((b-a)/dx)
x = np.linspace(a,b, N+1)
ind_space = (N*dxspace/(b-a)).astype(int)

#Discretize time
T = params.maxtime
maxD = max(D)
dt = 0.25*dx*dx/maxD
Nt = int(T/dt)
t = np.linspace(0, T, Nt+1)

F = dt/(dx*dx)

unew   = np.zeros(N+1) # unknown u at t+1
uold = np.zeros(N+1)   # u at t

# Set initial condition u(x,0) = I(x)
for i in range(0, N+1):
    u_1[i] = I(x[i])

for n in range(0, Nt):
    # Compute u at inner mesh points
    for i in range(1, Nx):
        u[i] = u_1[i] + F*(u_1[i-1] - 2*u_1[i] + u_1[i+1])

    # Insert boundary conditions
    u[0] = 0;  u[Nx] = 0

    # Update u_1 before next step
    u_1[:]= u