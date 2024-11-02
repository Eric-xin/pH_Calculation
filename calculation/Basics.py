import numpy as np

class Inert:
    def __init__(self, charge, conc=None):
        if charge is None:
            raise ValueError("The charge for this ion must be defined.")
        self.charge = charge
        self.conc = conc

    def alpha(self, pH):
        length = 1 if isinstance(pH, (int, float)) else len(pH)
        return np.ones(length).reshape(-1, 1)

# Actually 'Intert' is a subclass of 'Acid'. It is a special case of an acid with Ka = inf. Ka = inf means that the acid is completely dissociated in water. That is why the alpha of 'Inert' is always the same and equal to 1.


class Acid:
    def __init__(self, Ka=None, pKa=None, charge=None, conc=None):
        if Ka is None and pKa is None:
            raise ValueError("You must define either Ka or pKa values.")
        if charge is None:
            raise ValueError("The maximum charge for this acid must be defined.")

        if Ka is None:
            self.pKa = np.array([pKa] if isinstance(pKa, (int, float)) else pKa, dtype=float)
            self.pKa.sort()
            self.Ka = 10**(-self.pKa)
        else:
            self.Ka = np.array([Ka] if isinstance(Ka, (int, float)) else Ka, dtype=float)
            self.Ka.sort()
            self.Ka = self.Ka[::-1]
            self.pKa = -np.log10(self.Ka)

        self._Ka_temp = np.append(1., self.Ka)
        self.charge = np.arange(charge, charge - len(self.Ka) - 1, -1)
        self.conc = conc

    def alpha(self, pH):
        pH = np.array([pH] if isinstance(pH, (int, float)) else pH, dtype=float)
        h3o = 10.**(-pH)
        if len(h3o) > 1:
            h3o = np.repeat(h3o.reshape(-1, 1), len(self._Ka_temp), axis=1)

        power = np.arange(len(self._Ka_temp))
        h3o_pow = h3o**(power[::-1])
        Ka_prod = np.cumprod(self._Ka_temp)
        h3o_Ka = h3o_pow * Ka_prod

        den = h3o_Ka.sum(axis=1) if len(h3o.shape) > 1 else h3o_Ka.sum()
        return h3o_Ka / den.reshape(-1, 1) if len(h3o.shape) > 1 else h3o_Ka / den

# Implementation of the 'Base' class. 'Base' is a kind of 'Acid' where charge is positive. The difference is that 'Base' has a 'Kb' value instead of 'Ka'. The 'Kb' value is the base dissociation constant.

# def Base(Kb=None, pKb=None, charge=None, conc=None):
#     if Kb is None and pKb is None:
#         raise ValueError("You must define either Kb or pKb values.")
#     if charge is None:
#         raise ValueError("The maximum charge for this base must be defined.")
#     if pKb is not None and Kb is None:
#         if isinstance(pKb, (int, float)):
#             Kb = 10**(-pKb)
#         else:
#             Kb = 10**(-np.array(pKb, dtype=float))
    
#     if isinstance(Kb, (int, float)):
#         Ka = 1./Kb
#     else:
#         Ka = 1./np.array(Kb, dtype=float)
    
#     return Acid(Ka, -charge, conc)
    