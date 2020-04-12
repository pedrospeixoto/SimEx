#Micro extraction model 
# Pedro Peixoto (ppeixoto@usp.br)

#Load libraries
import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt
import sys
plt.rcParams.update({'font.size': 18})

#Main model file
import mextractmodel as mex

#Setup compartments based on mex_params.py
p=mex.device()

#Figure parameters
fig, axes = plt.subplots(1,1, constrained_layout=True, figsize=(15,10))
plt.xlabel("x-distance (cm)")
plt.ylabel("Concentration ($\mu g/mL$)")
plt.title("Microextration Model")
for i, name in enumerate(p.xnames):
    x=0.8*p.xspace[i]+0.2*p.xspace[i+1]
    plt.text(x, -0.1, name)

#plt.ylim(bottom=-0.15) 
#plt.ylim(bottom=-0.15) 

#loop over time
print("\n i    time  mass")
for i, t in enumerate(p.time):    

    #Plot if required
    #if i%iplot == 0:
    if t in p.iplot:
        print(i, t, p.mass)

        #Plot
        axes.plot(p.x, p.uext, label=str(t)+"s")
        axes.legend()
        istr="{:07.0f}".format(t)
        plt.savefig(p.basename+"_"+istr+".png")
        plt.pause(0.05)
        np.savetxt(p.basename+"_"+istr+".csv", p.uext, delimiter=',')

    #Run time step
    p.run_timestep(t=t)

axes.plot(p.x, p.u_equi_ext, 'k--', label="Theoretical \n Equilibrium")
axes.legend()
plt.savefig(p.basename+"_final.png")
plt.show()
