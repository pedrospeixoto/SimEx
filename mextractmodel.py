import numpy as np
import scipy.sparse as sparse
from scipy.sparse.linalg import spsolve
import sys
import mex_params as params

class device:
    def __init__(self):
        self.header='''
        --------------------------------------------------------------
        Micro Extraction Diffusion Model
        Pedro S. Peixoto - University of São Paulo (ppeixoto@usp.br)
        --------------------------------------------------------------
        '''
        print(self.header)
        #----------------------------------
        # Extraction mechanism parameters defined via param.py
        #----------------------------------
        self.D=params.D #diffusion coefficients
        self.K=params.K #partition coefficients
        self.C=params.C #initial concentrations
        self.xspace=params.x #Mechanism domain and interfaces points
        self.ncomp=len(self.D) #number of compartments
        self.nparts=len(self.K) #number of interfaces
        self.domain_len=self.xspace[-1]-self.xspace[0] #Domain size


        print("You defined a mechanism with "+str(self.ncomp)+" compartment(s).")
        print("Mechanism layout/interfaces (x): ",self.xspace)
        print("Initial concentrations:", self.C)
        print("Diffusion coefficients:", self.D)
        print("Interface coefficients:", self.K)
        print()

        #Check dimensions
        if self.nparts!=self.ncomp+1 and self.ncomp>1:
            print("Number of partitions must match the number of interfaces between spaces")
            print("Please re-configure parameters in params.py")
            sys.exit(-1)

        #Padd the D vectors with 0.0 at 1st and last positions
        #This is just to simplify knowledge of the boundaries 
        self.D = np.insert(self.D, 0, 0., axis=0)
        self.D = np.append(self.D,[0.])

        #Initialize compartment solvers
        #Create list of compartments
        self.compart = []
        for i in range(1, self.ncomp+1):
            #Initialize compartment
            Dloc=np.array([self.D[i-1],self.D[i], self.D[i+1]])
            Kloc=np.array([self.K[i-1],self.K[i]])
            xloc=np.array([self.xspace[i-1],self.xspace[i]])
            self.compart.append( self.compartment(i-1, Dloc, Kloc, xloc))

        #Discretize compartments and initialize the concentration on the grid
        self.Ninit=params.N
        print("Proposed number of control volumes (grid points): ", self.Ninit)
        self.N=0     
        for i, comp in enumerate(self.compart):
            ni=self.N
            n=int(self.Ninit*(comp.len/self.domain_len))-1
            comp.init_disc(n, ni) #configure 
            self.N=self.N+comp.n
            print("Compart: ", comp.icomp, " ini_index:", comp.ni, " deg_free",  comp.n)    
        self.ndf=self.N
        self.N=self.N+self.nparts #grid points, for plotting
        self.dx=(self.domain_len)/(self.N)
        print("Adjusted number of grid points: ", self.N)      
        print("Number of dregres of freedom: ", self.ndf)
        
        #Define global matrix
        main  = np.ones(self.ndf)
        lower = np.ones(self.ndf-1)
        upper = np.ones(self.ndf-1)

    
        #print(self.A.todense())
        for i, comp in enumerate(self.compart):
            comp.build_sys(main, lower, upper)
        #print(self.A.todense())
        #A=sparse.dia_matrix((main, 0), shape=(self.ndf, self.ndf))
        #A=A+sparse.dia_matrix((lower, -1), shape=(self.ndf, self.ndf))

        #Fill matrix with compartment information (pre-computation)
        self.A = sparse.diags(
            diagonals=[main, lower, upper],
            offsets=[0, -1, 1], shape=(self.ndf, self.ndf),
            format='csr')

        #print(self.A.todense())
        self.I=sparse.identity(self.ndf,  format='csr')

        #self.A = sparse.diags(
        #    diagonals=[main, lower, upper],
        #    offsets=[0, -1, 1], shape=(self.ndf, self.ndf),
        #    format='csr')

        #Fill initial conditions
        self.u = np.zeros(self.ndf)
        
        for i, comp in enumerate(self.compart):            
            self.u[comp.ni:comp.ni+comp.n]=np.full(comp.n, self.C[i])

        #print(self.u)
        self.extend_u()
        #print(self.uext)
        self.uold=np.copy(self.u)
        

    def extend_u(self):
        self.uext = np.copy(self.u)
        extramass=0
        for i, comp in enumerate(self.compart):
            #print(i, comp.n, comp.ni)
            if comp.K[0]==0:
                self.uext = np.insert(self.uext, comp.ni+i, self.u[comp.ni])
            else :
                uinter=(comp.D[1]/(comp.D[1]+comp.D[0]*comp.K[0]))*self.u[comp.ni]
                uinter=uinter+(comp.D[0]/(comp.D[1]+comp.D[0]*comp.K[0]))*self.u[comp.ni-1]
                extramass=extramass+uinter*comp.K[0]
                self.uext = np.insert(self.uext, comp.ni+i, uinter)
            #print(self.uext)
        self.mass=self.dx*(np.sum(self.uext)+extramass)
        return

    def precomp(self, dt):
        #Crank-nicolson matrices
        self.Bplus=self.I+(0.5*dt)*self.A
        self.Bminus=self.I-(0.5*dt)*self.A           
        #print(self.Bplus.todense())
        return

    def run_timestep(self, dt):
        #dt=dt
        b = self.Bplus.dot(self.u)
        #self.u = self.u+dt*self.A.dot(self.u) #spsolve(self.Bminus, b)
        self.u = spsolve(self.Bminus, b)
        #self.uold=self.u
        #print(self.Cold)
        #print(dt*self.A.todense())
        #print(dt*self.A.dot(self.Cold))
        #print(self.Cold+dt*self.A.dot(self.Cold))
        #sys.exit(1)
        return self.u #self.Cold+dt*self.A.dot(self.Cold)


    class compartment:

        def __init__(self, i,  D, K, x):
            self.icomp = i
            self.D=D
            self.K=K
            self.domain=x
            self.len=x[1]-x[0]
            print("Compartment setup")
            print("Local Domain:          ", x)
            print("Difusion (neigbours):  ", D)
            print("Border/Interfaces Coef:", K)
            print()
        
        def init_disc(self, n, ni):
            self.n=n #number of deg fredom
            self.ni=ni #start index in global matrix

        #Define space  domains
        def build_sys(self, main, lower, upper):
            #Space
            self.dx=(self.domain[1]-self.domain[0])/(self.n+1)
            #print(self.dx)

            # Precompute sparse matrix
            upper[self.ni:self.ni+self.n] = +(1/(self.dx*self.dx))*self.D[1]
            main[self.ni+1:self.ni+self.n-1]=-(2/(self.dx*self.dx))*self.D[1]
            lower[self.ni:self.ni+self.n-1] = +(1/(self.dx*self.dx))*self.D[1]

            #left border
            if self.K[0] != 0:
                main[self.ni] = -(1/(self.dx*self.dx))*self.D[1]*(self.D[1]+2.0*self.D[0]*self.K[0])/(self.D[1]+self.D[0]*self.K[0])
                lower[self.ni-1] = +(1/(self.dx*self.dx))*self.D[1]*(self.D[0])/(self.D[0]*self.K[0]+self.D[1])
            else:
                main[self.ni] = -(1/(self.dx*self.dx))*self.D[1]
            
            #right border
            if self.K[1] != 0:
                main[self.ni+self.n-1] = -(1/(self.dx*self.dx))*self.D[1]*(self.D[1]*self.K[1]+2.0*self.D[2])/(self.D[2]+self.D[1]*self.K[1])
                upper[self.ni+self.n-1] = +(1/(self.dx*self.dx))*self.D[1]*(self.D[2]*self.K[1])/(self.K[1]*self.D[1]+self.D[2])
            else:
                main[self.ni+self.n-1] = -(1/(self.dx*self.dx))*self.D[1]

            #print(main, lower, upper)            
            
            #np.copy(-self.dif[0:-1]-self.dif[1:]) #1
            #lower = np.copy(self.dif[1:-1]) #1
            #upper = np.copy(lower)
           
            
            #self.A=(1/self.dx*self.dx)*A
            #print(self.A.todense())    

            #Crankn-cholson matrices
            #self.I=sparse.identity(n,  format='csr')
            

        

    
        


    




