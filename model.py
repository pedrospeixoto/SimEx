#Micro extraction model 
# Pedro Peixoto (ppeixoto@usp.br)

import numpy as np 
import scipy as sp 
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys

import params
import compart_dif_solver 

import matplotlib.pyplot as plt


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
K=params.K #partition coefficients
C=params.C #initial concentrations
xspace=params.x #Mechanism domain and interfaces points
ncomp=len(D) #number of compartments
nparts=len(K) #number of interfaces
domain_len=xspace[-1]-xspace[0] #Domain size


print("You defined a mechanism with "+str(ncomp)+" compartment(s).")
print("Mechanism layout/interfaces (x): ",xspace)
print("Initial concentrations:", C)
print("Diffusion coefficients:", D)
print("Interface coefficients:", K)
print()

#Check dimensions
if nparts!=ncomp+1 and ncomp>1:
    print("Number of partitions must match the number of interfaces between spaces")
    print("Please re-configure parameters in params.py")
    sys.exit(-1)

#Padd the D vectors with 0.0 at 1st and last positions
#This is just to simplify knowledge of the boundaries 
D = np.insert(D, 0, 0., axis=0)
D = np.append(D,[0.])

#Initialize compartment solvers
#Create list of compartments
compart = []
for i in range(1, ncomp+1):
    #Initialize compartment
    Dloc=np.array([D[i-1],D[i], D[i+1]])
    Kloc=np.array([K[i-1],K[i]])
    xloc=np.array([xspace[i-1],xspace[i]])
    compart.append(compart_dif_solver.compartment(Dloc, Kloc, xloc))

#Discretize compartments and initialize the concentration on the grid
Ninit=params.N
print("Proposed number of space grid points: ", Ninit)
N=0
for i, comp in enumerate(compart):
    n=int(Ninit*(comp.len/domain_len))
    N=N+n
    comp.discretize(n)
    comp.init_conditions(C[i]) #If you need variable initial concentration, change this function
print("Adjusted number of gridpoints: ", N)

#Discretize time
T=params.maxtime
maxD = max(D)
dx = domain_len/N
dt = 100*0.25*dx*dx/maxD
Nt = int(T/dt)
time = np.linspace(0, T, Nt+1)

#Pre-compute crank-nic matrices
for comp in compart:
    comp.precomp(dt)

#Join init condition and full solution for plotting
C=compart_dif_solver.join_comparts_data(compart)

fig, axes = plt.subplots(1,1, constrained_layout=True, figsize=(10,15))
axes.plot(C)
plt.pause(0.05)


#loop over time
for i, t in enumerate(time):
    print(i,t)
    for comp in compart:
        comp.run_timestep(dt)
    C=compart_dif_solver.join_comparts_data(compart)
    axes.plot(C)
    plt.pause(0.05)

    
plt.savefig("microextraction_evolution.png")
plt.show()
