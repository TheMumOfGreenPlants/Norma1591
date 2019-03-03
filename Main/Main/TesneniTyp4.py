from math import *
from Tesneni import *
from Soucast import *

class TesneniTyp4(Tesneni):
    """description of class"""

    r_2 = 1         # polomer zakriveni v prurezu tesneni [mm]
    fi_G = 0        # uhel sklonu tesnici plochy [rad]

    def calcb_Gifirst(self): 
        """(74)"""
        self.b_Gi = (12 * self.r_2 * cos(self.fi_G) * self.b_Gt * self.Q_smax / self.E[0])**(1/2)

    def calcb_Gi(self,obj1,obj2,F_G0):
        """(75)"""
        self.b_Gi = ((12 * self.r_2 * cos(self.fi_G) * F_G0)/(pi * self.d_Ge * self.E_G0) + (F_G0 / (pi * self.d_Ge * self.Q_smax))**2)**(1/2)

    def calcd_Ge(self):
        """(76)"""
        self.d_Ge = self.d_Gt

    def calcE_Gm(self,F_G0):
        """prazdna funkce"""
        pass