import numpy as np
import scipy.sparse as sparse
from scipy.sparse.linalg import spsolve
import sys

def join_comparts_data(compart):
    C=np.array([])
    for comp in compart:
        C=np.concatenate((C, comp.C ))
    return C

class compartment:

    def __init__(self, D, K, x):
        self.D=D
        self.K=K
        self.domain=x
        self.len=x[1]-x[0]
        print("Compartment setup")
        print("Local Domain:          ", x)
        print("Difusion (neigbours):  ", D)
        print("Border/Interfaces Coef:", K)
        print()
        
    #Define space  domains
    def discretize(self, n=10):
        #Space
        self.n=n
        self.dx=(self.domain[1]-self.domain[0])/n
        self.x = np.linspace(self.domain[0],self.domain[1], n+1)
        self.xhalf = np.linspace(self.domain[0]+self.dx/2.0, self.domain[1]-self.dx/2.0, n)

        #Diffusion vector
        self.dif = np.full(self.n+1, self.D[1]) 
        #Override border with average
        self.dif[0]=0.5*(self.D[0]+self.D[1])
        self.dif[-1]=0.5*(self.D[1]+self.D[2])

        #Define problem matrix
        main  = np.zeros(n)
        lower = np.zeros(n-1)
        upper = np.zeros(n-1)

        # Precompute sparse matrix
        main = np.copy(-self.dif[0:-1]-self.dif[1:]) #1
        lower = np.copy(self.dif[1:-1]) #1
        upper = np.copy(lower)

        #Boundaries
        
        main[0]=main[0]+self.dif[0]*self.K[0]
        main[-1]=main[-1]+self.dif[-1]/self.K[1]
        
        A = sparse.diags(
            diagonals=[main, lower, upper],
            offsets=[0, -1, 1], shape=(n, n),
            format='csr')
        
        self.A=(1/self.dx*self.dx)*A
        print(self.A.todense())    
        #Crankn-cholson matrices
        self.I=sparse.identity(n,  format='csr')
        

    #Initialize with constant per compartment initial conditions
    def init_conditions(self, C):
        #self.Cinit=np.full(self.n, C)
        self.Cinit=np.random.rand(self.n)
        self.C=np.copy(self.Cinit)
        self.Cold=np.copy(self.Cinit)

    def precomp(self, dt):
        #Crank-nicolson matrices
        self.Bplus=self.I+(0.5*dt)*self.A
        self.Bminus=self.I-(0.5*dt)*self.A           
        #print(self.Bplus.todense())

    def run_timestep(self, dt):
        b = self.Bplus.dot(self.Cold)
        self.C = spsolve(self.Bminus, b)
        self.Cold=self.C
        #print(self.Cold)
        #print(dt*self.A.todense())
        #print(dt*self.A.dot(self.Cold))
        #print(self.Cold+dt*self.A.dot(self.Cold))
        #sys.exit(1)
        return self.C #self.Cold+dt*self.A.dot(self.Cold)
        


    




