from Priruba import *

class TocivaPrirubaSObrubou_Lemem(Priruba):
    """description of class"""

    e_L = 18       # vypocet dle 2 * A_L /(d_4 - d_6)

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

