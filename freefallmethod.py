import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
from scipy.optimize import curve_fit
from truncate import truncate

def f_stokes(Re):
    return 0.

def f_cliftgauvin(Re):
    return 0.2415*Re**0.687+Re/12.*0.42/(1+19018/(Re**1.16))

def f_flemmerbanks(Re):
    Ret=2.*Re
    return 10.**(0.261*Ret**0.369 - 0.105*Ret**0.431 - 0.124 /(1+(np.log(Ret)/np.log(10.))**2)) -1.

def f_cheng(Re):
    Ret=2.*Re
    return (1+0.27*Ret)**0.43+Ret/24.*0.47*(1-np.exp(-0.04*Ret**0.38)) - 1.

def wlogifuncfree(x,k,x0,gamma):
   f = (1+np.exp(-k*(x-x0)))**gamma
   return f / (1-f)**(1./4.)

def logifuncfree(x,k,x0,gamma):
   f = (1+np.exp(-k*(x-x0)))**gamma
   return f

class FreeFallMethod:
    def __init__(self,name, shortname, f):
        self.name = name  # Name of the method
        self.shortname=shortname
        self.f    = f     # f function so that C_d = 12/Re * f

    def delta(self,R):
        # Evaluate the delta function as a function of the virtual Reynolds number Re
        d=0.
        eps=0.000000001
        error=1.
        while(error>eps):
            dprev = d
            d = 1./(1+self.f(R*(1+d))) - 1.
            error = np.abs(d-dprev)
        return d

    def evaluate_chengwstar(self):
        self.cdcheng=np.zeros((self.n,))
        self.wstarcheng=np.zeros((self.n,))
        
        for i in range(self.n):
            self.cdcheng[i] = 432/self.dstar[i]**3 * (1+0.022*self.dstar[i]**3)**0.54 + 0.47*(1-np.exp(-0.15*self.dstar[i]**0.45))
            self.wstarcheng[i] = np.sqrt(4.*self.dstar[i]/(3.*self.cdcheng[i]))
            
    def evaluate(self,rarray):
        self.n=rarray.shape[0]
        
        self.rarray   = rarray
        self.rlog = np.log(self.rarray)                 # log(R) in an array

        self.deltavec = np.zeros((self.n,))
        for i in range(self.n):
            self.deltavec[i]=self.delta(rarray[i])
            
        self.reyarray = rarray*(1+self.deltavec)     # Re in an array
        self.mdelta   = -self.deltavec

        self.dstar = (36.*self.rarray)**(1./3.)     # Notation of Cheng 2009
        self.wstar = 2.* self.reyarray / self.dstar # dimensionless settling velocity
        
    def logistic_fit(self):

        # Fit as a function of virtual Reynolds
        self.wpopt, self.wpcov = curve_fit(wlogifuncfree, self.rlog, self.mdelta/(1-self.mdelta)**(1./4.), p0=[0.4, 0.5, -1.9]) # For weighted optimization
        print('    weighted opt parameters on R number', self.wpopt)
        print('    Parameters std on R number', self.wpcov)

        # Roundup parameters
        self.wpopt=[truncate(x,4) for x in self.wpopt ] 
        print('    Rounded opt parameters on R number', self.wpopt)
        self.fitdelta=np.zeros(self.n,)
        for i in range(self.n):
            self.fitdelta[i]=-logifuncfree(self.rlog[i],*self.wpopt)
        self.fiterror=(self.fitdelta - self.deltavec) / (1+self.deltavec) * 100.
        self.wstarfit = 2.* self.rarray * (1+self.fitdelta) / self.dstar # dimensionless settling velocity
            
        # Fit as a function of real Reynolds
        self.reypopt, self.reypcov = curve_fit(wlogifuncfree, np.log(self.reyarray), self.mdelta/(1-self.mdelta)**(1./4.), p0=[0.44, 1.64, -1.5]) # For weighted optimization
        print('    weighted opt parameters on Reynolds', self.reypopt)
        self.reyfitdelta=np.zeros(self.n,)
        for i in range(self.n):
            self.reyfitdelta[i]=-logifuncfree(np.log(self.reyarray[i]),*self.reypopt)
        self.reyfiterror=(self.reyfitdelta - self.deltavec) / (1+self.deltavec) *100.

    def logistic_delta(self,R):
        return -logifuncfree(np.log(R),*self.wpopt)
