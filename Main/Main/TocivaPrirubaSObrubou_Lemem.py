from Priruba import *
from IntegralniPriruba import *
from TesneniTyp1 import *

class TocivaPrirubaSObrubou_Lemem(Priruba):
    """description of class"""

    def __init__(self, krk_volba):
        self.krk = krk_volba
        self.skorepina = 1

    e_L = 18        # vypocet dle 2 * A_L /(d_4 - d_6)
    b_0 = 0         # sirka zkoseni (nebo zaobleni) tocive priruby
    d_6 = 1         # vnitrni prumer tocive priruby
    d_8 = 2         # vnejsi prumer lemu/obruby
    #POUZE priruba 12
    e_1 = 4
    e_2 = e_1



    def calc622(self):
        """(11)(12) *(13)* (14)(15) *(16)*"""
        self.calc6221()
        self.b_F = (self.d_8 - self.d_0) / 2
        self.d_F = (self.d_8 + self.d_0) / 2
        #self.e_F = 2 * self.A_F / (self.d_8 - self.d_0)
        self.b_L = (self.d_4 - self.d_6) / 2 - self.d_5e
        self.d_L = (self.d_4 + self.d_6) / 2
        #self.e_L = 2 * self.A_L / (self.d_4 - self.d_6)

    def calc623(self):
        self.e_P =self.e_F
        if self.krk == 1:
            IntegralniPriruba.calc6231(self)
        elif self.krk == 0:
            IntegralniPriruba.calc6232(self)

    def calc624(self):
        """(40)"""
        IntegralniPriruba.calc624(self)
        self.Z_L = 3 * self.d_L / (pi * self.b_L * self.e_L**3)

        

    def calcd_70(self):
        """(61)"""
        self.d_70 = min(max(objPrirubaX.d_7min,(self.d_Ge + self.chi * objPrirubaX.d_3e) / (1 + self.chi)), objPrirubaX.d_7max)

    def calcchi(self):
        """(62)"""
        chi = (objPrirubaX.Z_L * objPrirubaX.E[0]) / (objPrirubaX.Z_F * objPrirubaX.E[0])

    def calcd_7min(self):
        """(85)"""
        self.d_7min = self.d_6 + 2 * self.b_0

    def calcd_7max(self):
        """(86)"""
        self.d_7max = self.d_8

    def calch_GHL(self):
        """(87)(88)(89)"""
        self.h_G = (self.d_70 - self.objTesneni.d_Ge) / 2
        self.h_H = (self.d_70 - self.d_E) / 2
        self.h_L = (self.d_3e - self.d_70) / 2

    def calcW_L(self):
        """(150)"""
        self.W_L = (pi/2) * self.f_L * self.b_L * self.e_L**2

    def calcPhi_L(self):
        """(149)"""
        self.calcW_L()
        self.Phi_L = F_B * self.h_L / self.W_L

    def calcPhi_F(self,objTesneni):
        """(129) nebo (151)"""
        IntegralniPriruba.calcPhi_F(self)
        if isinstance(objTesneni,TesneniTyp1) and (self.objTesneni.d_G2 - self.d_7 > 0):
            self.Phi_F = (abs(F_Q + F_R) * self.h_H) / ((pi/4)* self.d_E * (self.f_E * min(self.e_E**2,self.e_F**2) + min(self.f_F * self.e_F**2, Qmax * (self.d_G2 - self.d_7**2) / 4)))
