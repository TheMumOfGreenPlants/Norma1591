from math import *
from Soucast import *
from Tesneni import *

class TesneniTyp2(Tesneni):
    """Kovova tesneni se zaoblenymi povrchy, jednoduchy dotyk - obrazky 3b),c)"""

    r_2 = 1             # polomer zakriveni v prurezu tesneni [mm]
    fi_G = 0      # uhel sklonu tesnici plochy [rad]
    d_G0 = 0

    def calcb_Gifirst(self): 
        """(69)"""
        self.calcd_Ge()
        self.b_Gi = (6 * self.r_2 * cos(self.fi_G) * self.b_Gt * self.Q_smax / self.E[0])**(1/2)

    def calcb_Gi(self, F_G0):
        """(70)"""
        self.b_Gi = ((6 * self.r_2 * cos(self.fi_G) * F_G0)/(pi * self.d_Ge * self.E_G[0]) + (F_G0 / (pi * self.d_Ge * self.Q_smax))**2)**(1/2)

    def calcd_Ge(self):
        """(71)"""
        self.d_Ge = self.d_G0
