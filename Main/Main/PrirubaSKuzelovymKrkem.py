from Priruba import *

class PrirubaSKuzelovymKrkem(Priruba):
    """description of class"""
    e_1 = 26.5     # nejmensi tloustka steny na tenkem konci krku   [mm]
    e_2 = 28.5     # tloustka steny na silnem konci krku            [mm]
    l_H = 40       # delka krku                                     [mm]
    Fi_S = 0       # natoceni pripojne skorepiny                    [rad]
    d_1 = 91.5     # stredni prumer krku na tenci strane            [mm]
    d_2 = 93.5     # stredni prumer krku na silnejsi strane         [mm]



    def calce_E(self):
        """(17)"""
        self.calcBeta()
        self.e_E = self.e_1 * (1 + ( self.Beta - 1 ) * self.l_H / ( ( self.Beta / 3 ) * sqrt ( self.d_1 * self.e_1 ) + self.l_H ) )

    def calce_D(self):
        """(18)"""
        self.calcBeta()
        self.e_D = self.e_1 * ( 1 + ( self.Beta - 1 ) * self.l_H / ( (self.Beta / 3)**4 * (self.d_1 * self.e_1)**2 + self.l_H**4 )**(1/4))
    
    def calcBeta(self):
        """(19)"""
        self.Beta = self.e_2 / self.e_1

    def calcd_E(self):
        """(20)"""
        self.d_E = ( min( ( self.d_1 - self.e_1 + self.e_E ) , ( self.d_2 + self.e_2 - self.e_E ) ) + max( ( self.d_1 + self.e_1 - self.e_E ) , ( self.d_2 - self.e_2 + self.e_E) ) ) / 2