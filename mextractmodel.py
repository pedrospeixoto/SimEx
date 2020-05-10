#-----------------------------------------
# Main microextraction model library
# P. Peixoto (ppeixoto@usp.br)
#----------------------------------------

#Libraries
import sys
import os

import numpy as np
import scipy.sparse as sparse
from scipy.sparse.linalg import spsolve
import pandas as pd

#Default parameters file 
import mex_params as params

#Main class for device information
class device:
    def __init__(self, \
        D = params.D, \
        K = params.K, \
        C = params.C, \
        xspace = params.x, \
        xnames = params.xnames, \
        name = params.name, \
        N = params.N, \
        dt = params.dt, \
        maxtime = params.maxtime, \
        iplot_time = params.iplot_time ):

        self.header='''
        --------------------------------------------------------------
        Micro Extraction Diffusion Model
        Pedro S. Peixoto - University of Sao Paulo (ppeixoto@usp.br)
        --------------------------------------------------------------
        '''
        print(self.header)
        #----------------------------------
        # Extraction mechanism parameters defined via mex_param.py
        #----------------------------------
        self.D = D #diffusion coefficients
        self.K = K #partition coefficients
        self.C = C #initial concentrations
        self.xspace = xspace #Mechanism domain and interfaces points
        self.xnames = xnames #Names of compartments

        self.ncomp=len(self.D) #number of compartments
        self.nparts=len(self.K) #number of interfaces
        self.domain_len=self.xspace[-1]-self.xspace[0] #Domain size

        #output directory settings
        self.dir="output"
        if not os.path.exists(self.dir):
            os.makedirs(self.dir)
        self.basename = name
        self.basedir=self.dir+"/"+self.basename
        if not os.path.exists(self.basedir):
            os.makedirs(self.basedir)
        
        self.basename=self.basedir+"/"+self.basename

        print("You defined a mechanism with "+str(self.ncomp)+" compartment(s).")
        print("Mechanism layout/interfaces (x): ",self.xspace)
        print("Initial concentrations:", self.C)
        print("Diffusion coefficients:", self.D)
        print("Interface coefficients:", self.K)
        print("Output basename:", self.basename)
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
            Dloc = np.array([self.D[i-1],self.D[i], self.D[i+1]])
            Kloc = np.array([self.K[i-1],self.K[i]])
            xloc = np.array([self.xspace[i-1],self.xspace[i]])
            self.compart.append( self.compartment(i-1, Dloc, Kloc, xloc, xnames[i-1]))

        #Discretize compartments and initialize the concentration on the grid
        self.Ninit = N
        print("Proposed number of control volumes (grid points): ", self.Ninit)
        self.N = 0     
        for i, comp in enumerate(self.compart):
            ni = self.N
            n = int(self.Ninit*(comp.len/self.domain_len))-1
            comp.init_disc(n, ni) #configure 
            self.N = self.N + comp.n
            #print("Compart: ", comp.icomp, " ini_index:", comp.ni, " deg_free",  comp.n)    
        self.ndf = self.N
        self.N = self.N+self.nparts #grid points, for plotting
        self.dx = (self.domain_len)/(self.N)
        
        self.x=np.linspace(self.xspace[0], self.xspace[-1], self.N, endpoint=True)
        #self.x=self.x[:-1]
        
        print("Adjusted number of grid points: ", self.N)      
        print("Number of degrees of freedom: ", self.ndf)
        
        #Define global tridiagonal matrix
        main  = np.ones(self.ndf)
        lower = np.ones(self.ndf-1)
        upper = np.ones(self.ndf-1)

        #print(self.A.todense())
        for i, comp in enumerate(self.compart):
            comp.build_sys(main, lower, upper)
        
        #Fill matrix with compartment information (pre-computation)
        self.A = sparse.diags(
            diagonals=[main, lower, upper],
            offsets=[0, -1, 1], shape=(self.ndf, self.ndf),
            format='csr')

        #print(self.A.todense())
        self.I=sparse.identity(self.ndf,  format='csr')

        #Fill initial conditions
        self.u = np.zeros(self.ndf)
        for i, comp in enumerate(self.compart):            
            self.u[comp.ni:comp.ni+comp.n]=np.full(comp.n, self.C[i])

        #Extend solution to interfaces and endpoints
        #self.extend_u()
        self.uext, self.mass = self.extend(self.u)
        self.mass_war = True
        self.mass_ini = self.mass

        #Time definition
        #Discretize time
        self.T = maxtime
        self.maxD = max(self.D)
        self.dt = dt #0.1 #0.1*dx/maxD #0.25*dx*dx/maxD
        self.Nt = int(self.T/self.dt)
        self.time = np.linspace(0, self.T, self.Nt+1)
        self.iplot = iplot_time
        print()
        print("Time-space info (dx, dt, Nt, maxD, dx/maxD):")
        print(self.dx, self.dt, self.Nt, self.maxD, self.dx/self.maxD)
        print()


        #Calculate equilibrium solution - reference
        self.equilibrium()
        self.diff_to_eq(0.0)
        
        #Precompute matrices
        self.Bplus = self.I+(0.5*self.dt)*self.A
        self.Bminus = self.I-(0.5*self.dt)*self.A    
        #self.B=self.I-(self.dt)*self.A    

        print("------------------------------------------------")
        print()

    def extend_u(self):

        #Add information on boundary points
        self.uext = np.copy(self.u)
        extramass=0
        for i, comp in enumerate(self.compart):
            #print(i, comp.n, comp.ni)
            if comp.K[0]==0:
                self.uext = np.insert(self.uext, comp.ni+i, self.u[comp.ni])
            else :
                uinter = (comp.D[1]/(comp.D[1]+comp.D[0]*comp.K[0]))*self.u[comp.ni]
                uinter = uinter + (comp.D[0]/(comp.D[1]+comp.D[0]*comp.K[0]))*self.u[comp.ni-1]
                extramass = extramass + uinter*comp.K[0]
                self.uext = np.insert(self.uext, comp.ni+i, uinter)
            #print(self.uext)
        #self.mass=self.dx*(np.sum(self.uext)+extramass)
        self.mass=self.dx*(np.sum(self.uext)+extramass)

        return self.uext

    def extend(self, u):

        #Add information on boundary points
        uext = np.copy(u)
        extramass=0
        for i, comp in enumerate(self.compart):
            #print(i, comp.n, comp.ni)
            if comp.K[0]==0:
                uext = np.insert(uext, comp.ni+i, u[comp.ni])
            else :
                uinter = (comp.D[1]/(comp.D[1]+comp.D[0]*comp.K[0]))*u[comp.ni]
                uinter = uinter + (comp.D[0]/(comp.D[1]+comp.D[0]*comp.K[0]))*u[comp.ni-1]
                extramass = extramass + uinter*comp.K[0]
                uext = np.insert(uext, comp.ni+i, uinter)
            #print(self.uext)
        #self.mass=self.dx*(np.sum(self.uext)+extramass)
        #add last point
        uext = np.append(uext, uext[-1])
        mass=self.dx*(np.sum(uext)+extramass)
        

        return uext, mass

    def diff_to_eq(self, time):
        #Calculates the diffence in the solution to the equilibrium

        self.eq_dif = self.u_equi - self.u
        try:
            self.eq_perc = self.eq_dif/self.u_equi
        except:
            print("Equilitrium has null values, can't calulate percent error")
            sys.exit(1)

        self.eq_dif_max_abs = np.max(np.abs(self.eq_dif))
        self.eq_perc_max = np.max(np.abs(self.eq_perc))

        #Collect time when it reaches a certain percent of equilibrium
        for i in range(len(self.equi_percents)):
            if self.eq_perc_max < np.abs(1-self.equi_percents[i]):
                if time < self.equi_percents_times[i]:
                    self.equi_percents_times[i]=time
                    #print(self.eq_dif_max_abs, self.eq_perc_max, self.equi_percents, self.equi_percents_times)

        #Check if mass is not changing much
        delta_mass = np.abs(self.mass - self.mass_ini )/self.mass_ini
        if delta_mass > 0.01 and self.mass_war:
            print("Warning: mass changes are large, it may be better to increase N and reduce dt ", self.mass, self.mass_ini, delta_mass)
            self.mass_war = False #only give 1 warning

    def run_timestep(self, t=0.0):
        #self.u = self.u+dt*self.A.dot(self.u) #Euler scheme
        #self.u = spsolve(self.B, self.u) #Implicit Euler
        self.u = spsolve(self.Bminus, self.Bplus.dot(self.u)) #Crank-Nicolson

        #self.extend_u()
        self.uext, self.mass = self.extend(self.u)
        self.diff_to_eq(t)
        
        return self.u 

    def run(self):
        
        u_snapshots = []
        u_df = pd.DataFrame()
        u_df["x"] = self.x
        for i, t in enumerate(self.time):    
            #Save when required
            #if i%iplot == 0:
            if t in self.iplot:
                
                mass_str = "{:.4f}".format(self.mass)
                perc_str = "{:.4f}".format(self.eq_perc_max*100.0)
                print(" It: ", i, " Time: ", t, " Mass: ", mass_str , " %Dif Eq: ", perc_str, "%" )

                #Plot
                u_snapshots.append(self.uext)                 
                istr="{:07.0f}".format(t)
                u_df["t"+istr] = self.uext
                filename = self.basename+"_"+istr+".csv"
                np.savetxt(filename, self.uext, delimiter=',')
                #print("  Saving snapshot at time "+istr+" with name: \n  ", filename)
                #print(np.average(self.uext))

            #Run time step
            self.run_timestep(t=t)

        #Save full matrix to file
        filename = self.basename+"_data.csv"
        print(" Saving data file with all times as ", filename)
        u_df.to_csv(filename)
        self.print_output()

        return u_snapshots

    def print_output(self):
        print()
        print(" %Eq   Time " )
        for i in range(len(self.equi_percents)):
            if self.equi_percents_times[i] > self.T : 
                time = "Not reached"
                print(self.equi_percents[i]*100,"%  ", time)
            else:
                time = "{:.2f}".format(self.equi_percents_times[i])
                print(self.equi_percents[i]*100,"%  ",  time)
        print()

    def equilibrium(self):
        #Calculate equilibrium solution

        #Initial mass
        Mini=0.0
        K=np.copy(self.K)
        K[0]=1.0
        for i, c in enumerate(self.C):
            dx=self.xspace[i+1]-self.xspace[i]
            Mini=Mini+c*dx
        self.Mini=Mini
        #print("Initial Mass:", Mini, np.sum(self.uext)*self.dx, self.mass)

        #Final distribution on compart 0
        a=0.0
        for i in range(self.ncomp):
            dx=self.xspace[i+1]-self.xspace[i]
            Kprod=1.0
            for k in K[0:i+1]:
                Kprod=Kprod*k
            a = a + dx/Kprod
            #print(i, dx, Kprod, a)
        C0=Mini/a

        Cend = [C0]
        for i in range(self.ncomp-1):
            #print(i, K[i], K[i+1], Cend)
            Cend.append(Cend[i]/K[i+1])
        print("Equilibrium concentrations:\n", Cend)
        self.Cend=Cend
        Mend=0.0
        for i, c in enumerate(Cend):
            dx=self.xspace[i+1]-self.xspace[i]
            Mend=Mend+c*dx

        self.u_equi = np.zeros(self.ndf)
        for i, comp in enumerate(self.compart):            
            self.u_equi[comp.ni:comp.ni+comp.n]=np.full(comp.n, self.Cend[i])
        
        self.u_equi_ext, self.mass_equi = self.extend(self.u_equi)
        
        #Control solution with respect to equilibrium
        self.equi_percents = [ 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.995] 
        self.equi_percents_times = [9999999999.0]*len(self.equi_percents)

    class compartment:

        def __init__(self, i,  D, K, x, name=""):
            self.icomp = i
            self.D=D
            self.K=K
            self.domain=x
            self.len=x[1]-x[0]
            print("Compartment", i, " setup")
            print(" Name:                  ", name)
            print(" Local Domain:          ", x)
            print(" Difusion (neigbours):  ", D)
            print(" Border/Interfaces Coef:", K)
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

            

        

    
        


    




