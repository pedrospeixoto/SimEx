import numpy as np 

#Micro extraction mechanism parameters

# All parameters must be adopted in the same units of space/time/mass/volume
# here we adopt:
# - space units: "cm" 
# - time units: "seconds" (s)
# - mass: micro-grams (mug)
# - volume: cm^3 (mL)

#Domain definition (position of interfaces) 
# x=np.array([initial_position, interface_1, interface_2, final_position])
#x=np.array([0, 1]) # cm
x=np.array([0, 50, 100]) # cm

#Initial concentrations in each space
#  Assumed constant for each space
C=np.array([1.0, 0.0]) #mug/mL

#Diffusion coefficients for each space - manually defined
D=np.array([1.0, 2.0]) # 2 spaces (cm^2/s)

#Partition coefficients K for each interface - manually defined 
K=np.array([1]) # 1 interface, non-dimensional




#Max time definition
maxtime = 10*60 #10 minutes