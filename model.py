#Micro extraction model 
# Pedro Peixoto (ppeixoto@usp.br)

#Load libraries
import numpy as np 
import scipy as sp 
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
plt.rcParams.update({'font.size': 18})

#Main model file
import mextractmodel as mex

#Main parameters files
import mex_params as params

#Setup compartments based on mex_params.py
p=mex.device()

#Discretize time
T=params.maxtime
maxD = max(p.D)
dx = p.domain_len/p.N
dt = 0.1 #0.1*dx/maxD #0.25*dx*dx/maxD
Nt = int(T/dt)
time = np.linspace(0, T, Nt+1)
print("\nTime-space info (dx, dt, Nt, maxD, dx/maxD):\n", dx, dt, Nt, maxD, dx/maxD, "\n")
iplot=params.iplot_time

#Figure parameters
fig, axes = plt.subplots(1,1, constrained_layout=True, figsize=(15,10))
plt.xlabel("x-distance (cm)")
plt.ylabel("Concentration ($\mu g/mL$)")
plt.title("Microextration Model")
for i, name in enumerate(params.xnames):
    x=0.8*params.x[i]+0.2*params.x[i+1]
    plt.text(x, -0.1, name)

#plt.pause(0.05)
#plt.show()


#loop over time
print("i    time     mass")
for i, t in enumerate(time):    

    #Plot if required
    #if i%iplot == 0:
    if t in iplot:
        print(i,t, p.mass)

        #Plot
        axes.plot(p.x, p.uext, label=str(t)+"s")
        axes.legend()
        istr="{:07.0f}".format(t)
        plt.savefig(p.basename+"_"+istr+".png")
        #plt.pause(0.05)

    #Run time step
    p.run_timestep(dt)


#plt.show()
