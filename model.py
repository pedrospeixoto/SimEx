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
nspaces=len(D) #number of spaces
nparts=len(K)
print("You defined a mechanism with "+str(nspaces)+" space(s).")
print("Diffusion coefficients:", D)
print("Partition coefficients:", K)

#Check dimensions
if nparts!=nspaces-1 and nspaces>1:
    print("Number of partitions must match the number of interfaces between spaces")
    print("Please re-configure parameters in params.py")
    sys.exit(-1)



