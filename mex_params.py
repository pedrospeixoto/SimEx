#Micro extraction mechanism parameters
import numpy as np 

# All parameters must be adopted in the same units of space/time/mass/volume
# here we adopt:
# - space units: "cm" 
# - time units: "seconds" (s)
# - mass: micro-grams (mug)
# - volume: cm^3 (mL)

#Domain definition (position of interfaces) 
# x=np.array([initial_position, interface_1, interface_2, final_position])
#x=np.array([0, 1]) # cm
# x0   x1   x2
# |----|----|
x=np.array([0, 0.5, 1]) # cm

#Initial concentrations in each space
#  Assumed constant for each space
#   C0   C1
# |----|----|
C=np.array([0.0, 1.0]) #mug/mL

#Diffusion coefficients for each space - manually defined
#   D0   D1
# |----|----|
D=np.array([2.0, 1.0]) # 2 spaces (cm^2/s)

#Partition coefficients K for each interface - manually defined 
# Set 1.0 for beggining and endpoints    
# K0  K1  K2
# |---|---|
K=np.array([0.0, 2.0, 0.0]) # 2 boundaries and 1 interface coefficient, non-dimensional

#Max time definition
maxtime = 100*60 #10 minutes

#Plotting time interval (every XX time units)
iplot_time=10

#Space discretization (number of grid points)
N = 10

#Experiment basename
name = "microextraction"
