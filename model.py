#Micro extraction model 
# Pedro Peixoto (ppeixoto@usp.br)

import numpy as np 
import scipy as sp 
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys

#Main model file
import compart_dif_solver 
#Main parameters files
import params

#Setup compartments based on params.py
p=compart_dif_solver.compartments()

#Discretize time
T=params.maxtime
maxD = max(p.D)
dx = p.domain_len/p.N
dt = 100*0.25*dx*dx/maxD
Nt = int(T/dt)
time = np.linspace(0, T, Nt+1)

#Pre-compute crank-nic matrices
for comp in p.compart:
    comp.precomp(dt)

#Join init condition and full solution for plotting
C=compart_dif_solver.join_comparts_data(p.compart)

fig, axes = plt.subplots(1,1, constrained_layout=True, figsize=(10,15))
axes.plot(C)
plt.pause(0.05)


#loop over time
for i, t in enumerate(time):
    print(i,t)
    for comp in p.compart:
        comp.run_timestep(dt)
    C=compart_dif_solver.join_comparts_data(p.compart)
    axes.plot(C)
    plt.pause(0.05)

    
plt.savefig("microextraction_evolution.png")
plt.show()
