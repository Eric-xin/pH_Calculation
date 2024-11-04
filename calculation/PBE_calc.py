import numpy as np
from scipy.optimize import minimize
# from PBE_Acids import *

try:
    from .PBE_Acids import *
except ImportError:
    from PBE_Acids import *

class PBE_calc:
    def __init__(self, *species, Kw=1.01e-14):
        self.species = species
        self.Kw = Kw
    
    def PBE_error(self, pH):
        pH = np.atleast_1d(pH).astype(float)
        h3o = 10.**(-pH)
        oh = self.Kw / h3o
        P_error = h3o - oh # H2O/H+ (Gain) H2O/OH- (Loss)

        for s in self.species:
            P_error += (s.conc * s.proton * s.alpha(pH)).sum(axis=-1)
        
        return np.abs(P_error)
    
    def pH_calc(self, guess=7.0, guess_est=False, est_num=1500, method='Nelder-Mead', tol=1e-5):
        if guess_est:
            phs = np.linspace(0, 14, est_num)
            guess = phs[self.PBE_error(phs).argmin()]
            
        self.pH_eq = minimize(self.PBE_error, guess, method=method, tol=tol)
        
        if not self.pH_eq.success:
            Warning('pH calculation did not converge')
            print(self.pH_eq.message)
            
        self.pH = self.pH_eq.x[0] if self.pH_eq.x.size == 1 else self.pH_eq.x