from Priruba import *

class TocivaPrirubaSObrubou_Lemem(Priruba):
    """description of class"""

    e_L = 18        # vypocet dle 2 * A_L /(d_4 - d_6)
    b_0 = 0         # sirka zkoseni (nebo zaobleni) tocive priruby
    d_6 = 1         # vnitrni prumer tocive priruby
    d_8 = 2         # vnejsi prumer lemu/obruby

    def calcbde_FL(self):
        """(11)(12) *(13)* (14)(15) *(16)*"""
        self.b_F = (self.d_8 - self.d_0) / 2
        self.d_F = (self.d_8 + self.d_0) / 2
        self.b_L = (self.d_4 - self.d_6) / 2 - self.d_5e
        self.d_L = (self.d_4 + self.d_6) / 2

    def calcfromIP(self):
        """(25)-(35)"""
        IntegralniPriruba.calcGama(self)
        IntegralniPriruba.calcTheta(self)
        IntegralniPriruba.calcLambda(self)
        IntegralniPriruba.calcGama(self)

    def calcZ_L(self):
        """(40)"""
        self.calcfromIP()
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
        self.h_G = (self.d_70 - Tesneni.d_Ge) / 2
        self.h_H = (self.d_70 - self.d_E) / 2
        self.h_L = (self.d_3e - self.d_70) / 2
