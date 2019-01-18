from Priruba import *

class ZaslepovaciPriruba(Priruba):
    """description of class"""

    def calcbde_FL(self):
        IntegralniPriruba.calcbde_FL()

    def calce_E(self):
        """(23)""" 
        self.e_E = 0

    def calcd_E(self):
        """(24)"""
        self.d_E = self.d_0
    
    def calcro(self):
        """(36)"""
        self.ro = self.d_9 / self.d_E

    def calch_R(self):
        """(37)"""
        self.calcro()
        self.h_R = (self.d_E / 4) * (1 - self.ro**2) * (0.7 + 3.3 * self.ro**2) / ((0.7 + 1.3 * self.ro**2) * (1 + self.ro**2))

    def calcZ_FL(self):
        """(38) (39)"""
        self.Z_F = 3 * self.d_F / (pi * (self.b_F * self.e_F**3 + self.d_F * self.e_0**3 * (1 - self.ro**2) / (1.4 + 2.6 * self.ro**2)))
        self.Z_L = 0

    def calce_P(self):
        """(78)"""
        self.e_P = 0