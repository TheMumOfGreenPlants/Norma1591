#from math import pi, sqrt, cos
from math import *
from Soucast import *
from Tesneni import *


class Priruba(Soucast):
    """description of class zkouska"""

# VSTUPNI PARAMETRY
    d_0 = 400
    d_3 = 495
    d_4 = 540
    d_5 = 23
    d_5t = 0
    l_5t = 0
    e_Fb = 36
    e_S = 0
    e_Ft = 39.5
    e_F = 36        # vypocet dle 2 * A_F /(d_4 - d_0)
    Fi_S = 0        # 0 - pro valec, natoceni pripojne skorepiny                    [rad]
    #POUZE priruba bez krku
    d_S = 403        # stredni prumer skorepiny (prumer v miste spoje s prirubou)        [mm]
    #E_F0 = 200000
    f_F = 205         # dovolene namahani priruby     [MPa]
    f_S = 205         # dovolene namahani skorepiny   [MPa]
# KONEC - VSTUPNI PARAMETRY

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

    def calch_P(self, d_Ge):
        """(77)"""
        self.h_P = ((d_Ge - self.d_E)**2 * (2 * d_Ge + self.d_E) / 6 + 2 * self.e_P**2 * self.d_F) / d_Ge**2

    def VypocitejPrirubu(self):
        # promenne totozne pro obe priruby
        self.beforecalc()
        self.calc622()
        self.calc623()
        self.calc624()