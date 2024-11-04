import numpy as np
from scipy.optimize import minimize
# from .CBE_Acids import *
# from CBE_Acids import *

try:
    from .CBE_Acids import *
except ImportError:
    from CBE_Acids import *

class CBE_calc:
    def __init__(self, *species, Kw=1.01e-14):
        self.species = species
        self.Kw = Kw
    
    def Charge_diff(self, pH):
        pH = np.atleast_1d(pH).astype(float) # Ensure pH is an array
        h3o = 10.**(-pH)
        oh = self.Kw / h3o
        x = h3o - oh

        for s in self.species:
            x += (s.conc * s.charge * s.alpha(pH)).sum(axis=-1)
        
        return np.abs(x)
    
    def pH_calc(self, guess=7.0, guess_est=False, est_num=1500, method='Nelder-Mead', tol=1e-5):
        if guess_est:
            phs = np.linspace(0, 14, est_num)
            guess = phs[self.Charge_diff(phs).argmin()]
            
        self.pH_eq = minimize(self.Charge_diff, guess, method=method, tol=tol)
        
        if not self.pH_eq.success:
            Warning('pH calculation did not converge')
            print(self.pH_eq.message)
            
        self.pH = self.pH_eq.x[0] if self.pH_eq.x.size == 1 else self.pH_eq.x
