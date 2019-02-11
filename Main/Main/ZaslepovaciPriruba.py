from Priruba import *
from IntegralniPriruba import *

class ZaslepovaciPriruba(Priruba):
    """description of class"""

    e_0 = 10
    d_9 = 0

    def calc622(self):
        IntegralniPriruba.calc622(self)

    def calc623(self):
        """(23)(24)""" 
        self.e_E = 0
        self.d_E = self.d_0

    def calc624(self):
        """(36)(37)(38)(39)"""
        self.ro = self.d_9 / self.d_E
        self.h_R = (self.d_E / 4) * (1 - self.ro**2) * (0.7 + 3.3 * self.ro**2) / ((0.7 + 1.3 * self.ro**2) * (1 + self.ro**2))
        self.Z_F = 3 * self.d_F / (pi * (self.b_F * self.e_F**3 + self.d_F * self.e_0**3 * (1 - self.ro**2) / (1.4 + 2.6 * self.ro**2)))
        self.Z_L = 0

    def calce_P(self):
        """(78)"""
        self.e_P = 0

    def calch_Q(self):
        """(80)"""
        self.h_Q = (self.d_E / 8) * (1 - self.ro**2) * (0.7 + 3.3 * self.ro**2) / (0.7 + 1.3 * self.ro**2) * (self.d_E / self.objTesneni.d_Ge)**2

    def calcPhi_F(self):
        """(145)"""
        self.calcW_F()
        self.Phi_F = max(abs(F_B * self.h_G + F_Q * (1 - self.ro**3) * self.objTesneni.d_Ge / 6 + F_R * (1 - self.ro) * self.objTesneni.d_Ge / 2),\
           abs(F_B * self.h_G + F_Q * (1 - self.ro**3) * self.objTesneni.d_Ge / 6), abs(F_R * (1 - self.ro) * self.objTesneni.d_Ge / 2)) / self.W_F

    def calcW_F(self):
        """(146)"""
        self.W_F = (pi/4) * self.f_F * (2 * self.b_F * self.e**2 + self.d_0 * (1 - self.ro) * self.e_0**2)

    def calcPhi_X(self):
        """(147)"""
        self.calcW_X
        self.Phi_X = F_B * (self.d_3 - self.d_X) / (2 * self.W_X)

    def calcW_X(self):
        """(148)"""
        self.W_X = (pi/4) * self.f_F * ((self.d_4 - 2 * self.d_5e - self.d_X) * self.e**2 + self.d_X * self.e_X**2)

    
