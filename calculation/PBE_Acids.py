import numpy as np

class PBE_Inert:
    def __init__(self, proton=None, proton_ref=None, conc=None):
        # if proton is None:
        #     raise ValueError("The proton for this ion must be defined.")
        if proton_ref is None:
            raise ValueError("The reference proton for PBE must be defined.")
        # self.proton = proton - proton_ref
        self.proton = -proton_ref
        self.conc = conc

    def alpha(self, pH):
        length = 1 if isinstance(pH, (int, float)) else len(pH)
        return np.ones(length).reshape(-1, 1)

class PBE_Acid:
    def __init__(self, Ka=None, pKa=None, proton=None, proton_ref=None, conc=None):
        if Ka is None and pKa is None:
            raise ValueError("You must define either Ka or pKa values.")
        if proton is None:
            raise ValueError("The maximum proton for this acid must be defined.")
        if proton_ref is None:
            raise ValueError("The reference proton for PBE must be defined.")

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
        self.proton = np.arange(proton, proton - len(self.Ka) - 1, -1)
        self.proton = self.proton - proton_ref
        self.proton_ref = proton_ref
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