from Priruba import *
from IntegralniPriruba import *
from TesneniTyp1 import *

class TocivaPrirubaSObrubou_Lemem(Priruba):
    """description of class"""

    E_L = numpy.asarray([200000,200000])

    def __init__(self, krk_volba):
        self.krk = krk_volba
        self.skorepina = 1

    e_L = 18        # vypocet dle 2 * A_L /(d_4 - d_6)
    b_0 = 0         # sirka zkoseni (nebo zaobleni) tocive priruby
    d_6 = 140         # vnitrni prumer tocive priruby
    d_8 = 180        # vnejsi prumer lemu/obruby
    #POUZE priruba 12
    e_1 = 4
    e_2 = e_1
    l_H = 40
    d_1 = 91



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
        A = {
            1 : IntegralniPriruba.calc6231,
            0 : IntegralniPriruba.calc6232,
            }[self.krk]
        A(self)


    def calc624(self):
        """(40)"""
        IntegralniPriruba.calc624(self)
        self.Z_L = 3 * self.d_L / (pi * self.b_L * self.e_L**3)

    def calc645(self,d_Ge):
        self.calch_P(d_Ge)
        IntegralniPriruba.calch_Q(self,d_Ge)

        ##zkontrolovat, zda mam vsechny veliciny
        self.calch_GHL(d_Ge)


    def calch_G0(self,d_Ge):
        """(60)"""
        self.calcd_70(d_Ge)
        self.h_G0 = (self.d_70 - d_Ge) / 2

    def calcd_70(self,d_Ge):
        """(61)"""
        self.calcd_7min()
        self.calcd_7max()
        self.calcchi()
        self.d_70 = numpy.minimum(numpy.maximum(self.d_7min,(d_Ge + self.chi * self.d_3e) / (1 + self.chi)), self.d_7max)

    def calcchi(self):
        """(62)"""
        self.chi = (self.Z_L * self.E[0]) / (self.Z_F * self.E_L[0])

    def calcd_7min(self):
        """(85)"""
        self.d_7min = self.d_6 + 2 * self.b_0

    def calcd_7max(self):
        """(86)"""
        self.d_7max = self.d_8

    def calch_GHL(self,d_Ge):
        """(87)(88)(89)"""
        self.h_G = (self.d_70 - d_Ge) / 2
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
