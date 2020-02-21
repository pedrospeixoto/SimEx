#Micro extraction model 
# Pedro Peixoto (ppeixoto@usp.br)

import numpy as np 
import scipy as sp 
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys

#Main model file
import mextractmodel as mex

#Main parameters files
import mex_params as params

#Setup compartments based on params.py
p=mex.device()

#Discretize time
T=params.maxtime
maxD = max(p.D)
dx = p.domain_len/p.N
dt = 0.1*dx/maxD #0.25*dx*dx/maxD
Nt = int(T/dt)
time = np.linspace(0, T, Nt+1)

#Pre-compute crank-nic matrices
p.precomp(dt)

fig, axes = plt.subplots(1,1, constrained_layout=True, figsize=(10,15))
axes.plot(p.uext)
plt.pause(0.05)
#print(p.u)
#print(p.uext)

#loop over time
for i, t in enumerate(time):
    p.run_timestep(dt)
    if i%10 == 0:
        p.extend_u()
        print(i,t, p.mass) #, p.u, p.uext)
        axes.plot(p.uext)
        plt.pause(0.05)

    
plt.savefig("microextraction_evolution.png")
plt.show()
