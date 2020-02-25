#Micro extraction model 
# Pedro Peixoto (ppeixoto@usp.br)

#Load libraries
import numpy as np 
import scipy as sp 
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
mpl.use('Agg')

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
dt = 0.1*dx/maxD #0.25*dx*dx/maxD
Nt = int(T/dt)
time = np.linspace(0, T, Nt+1)
iplot=params.iplot_time

fig, axes = plt.subplots(1,1, constrained_layout=True, figsize=(15,10))
plt.xlabel("Device position (cm)")
plt.ylabel("Concentration")

#plt.pause(0.05)
#plt.show()
print("i    time     mass")

#loop over time
for i, t in enumerate(time):    
    if i%iplot == 0:
        print(i,t, p.mass) #, p.u, p.uext)
        axes.plot(p.x, p.uext)
        plt.title("Microextration t="+str(t) )
        #plt.pause(0.05)
        istr="{:07.0f}".format(i)
        plt.savefig(p.basename+"_"+istr+".png")
    p.run_timestep(dt)
#plt.show()
plt.savefig(p.basename+"_"+str(t)+".png")