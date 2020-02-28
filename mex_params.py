#Micro extraction mechanism parameters
import numpy as np 

# All parameters must be adopted in the same units of space/time/mass/volume
# here we adopt:
# - space units: "cm" 
# - time units: "seconds" (s)
# - mass: micro-grams (mug)
# - volume: cm^3 (mL)

#Experiment basename
#name = "test1"
name = "microextraction_zhang1993_fig2"


#Domain definition (position of interfaces) 
# x=np.array([initial_position, interface_1, interface_2, final_position])
#x=np.array([0, 1]) # cm
# x0   x1   x2
# |----|----|
#x=np.array([0, 0.5, 1]) # cm
x=np.array([0, 0.0056, 0.07, 0.1]) # cm Zhang 1993 Fig 2

#Names of compartments
#xnames = ["", "", ""]
xnames = ["coating", "headspace", "aqueous"]

#Initial concentrations in each space
#  Assumed constant for each space
#   C0   C1   C2
# |----|----|----|
#C=np.array([0.0, 1.0]) #mug/mL
C=np.array([0.0, 0.0, 1.0]) #mug/mL Zhang 1993 Fig 2

#Diffusion coefficients for each space - manually defined
#   D0   D1   D2
# |----|----|----|
#D=np.array([2.0, 1.0]) # 2 spaces (cm^2/s)
D=np.array([2.8e-6, 0.077, 1.8e-5]) # (cm^2/s) Zhang 1993 Fig 2

#Partition coefficients K for each interface - manually defined 
# Set 0.0 for beggining and endpoints    
# K0  K1  K2   K3
# |---|---|----|
#K=np.array([0.0, 2.0, 0.0]) # 2 boundaries and 1 interface coefficient, non-dimensional
K=np.array([0.0, 50.0, 0.2, 0.0]) # Zhang 1993 Fig 2

#Max time definition 
#maxtime = 10*60 #10 minutes
maxtime = 60 #seconds

#Time step size
dt = 0.01  #Seconds 

#Plotting time spots
iplot_time=np.array([0.0, 3.0, 15.0, 30.0, 60.0, 120])
#iplot_time=np.linspace(0, 50, 21, endpoint=True)

#Space discretization (number of grid points)
N = 500
