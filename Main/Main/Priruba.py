#from math import pi, sqrt, cos
from math import *
from Soucast import *
from Tesneni import *


class Priruba(Soucast):
    """description of class zkouska"""
# vstupy
    d_0 = 65
    d_3 = 170
    d_4 = 220
    d_5 = 26
    d_5t = 24
    l_5t = 31
    e_Fb = 31
    e_S = 26.5
    e_Ft = 34
    e_F = 31       # vypocet dle 2 * A_F /(d_4 - d_0)

    #E_F0 = 200000
    f_F = 1         # dovolene namahani priruby     [MPa]
    f_S = 1         # dovolene namahani skorepiny   [MPa]

# vypocty
    def sete(self):
        self.e = self.e_Ft

    def setn_B(self, valn_B):
        self.n_B = valn_B

    def calc6221(self):
        """(3)(4)(6)"""
        self.p_B = pi*self.d_3 / self.n_B
        self.d_5e = self.d_5 * sqrt(self.d_5/self.p_B)
        self.d_3e = self.d_3 * (1 - 2 / self.n_B**2)

    # def calcd_5(self):
    #    """(5)"""
    #    self.d_5 = d_5t * l_5t/e_Fb

    def calch_P(self, tesneni):
        """(77)"""
        self.h_P = ((tesneni.d_Ge - self.d_E)**2 * (2 * tesneni.d_Ge + self.d_E) / 6 + 2 * self.e_P**2 * self.d_F) / tesneni.d_Ge**2