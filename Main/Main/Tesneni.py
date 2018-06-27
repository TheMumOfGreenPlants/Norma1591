from math import *
from PrirubaSKuzelovymKrkem import *

class Tesneni(object):
    """description of class"""
    Q_A = 200       # priloha G - neni pozadovana mira netesnosti   [MPa]
    d_G1 = 67       # teoreticky vnitrni prumer tesnici plochy      [mm]
    d_G2 = 120      # teoreticky vnejsi prumer tesnici plochy       [mm]
    e_G = 2         # tloustka tesneni v nezatizenem stavu          [mm]
    Q_smax = 480    # maximalni dovoleny tlak na tesneni            [MPa]

    ##vypoctove parametry - uzivatel nemeni

    def calcb_Gifirst(self):        
        self.calcb_Gt()
        self.b_Gi = self.b_Gt
        

    def calcb_Gt(self):
        """(51)"""
        self.b_Gt = (self.d_G2 - self.d_G1)/2

    def geth_G0(self, objPrirubaX):
        """(59)"""
        h_G0 = ( objPrirubaX.d_3e - self.d_Ge ) / 2
        return h_G0

    def getb_Gi(self):
        self.calcb_Gi()
        return self.b_Gi

    def calcb_Gi(self):
        """(65)"""
        self.calcd_Ge()
        self.calcE_Gm()
        #self.b_Gi = self.b_Gt
        self.b_Gi = sqrt( self.e_G / ( pi * self.d_Ge * self.E_Gm) / ( self.geth_G0(self.objPrvniPriruba) * self.objPrvniPriruba.Z_F / self.objPrvniPriruba.E_F0
                                                            + self.geth_G0(self.objDruhaPriruba) * self.objDruhaPriruba.Z_F / self.objDruhaPriruba.E_F0 ) + 
                    ( self.F_G0 / pi * self.d_Ge * self.Q_smax))
        
    def calcb_Ge(self):
        """(55)"""
        self.calcb_Gt()
        self.calcb_Gifirst()
        #self.calcA_Ge()
        
        self.b_Ge = min( self.b_Gi, self.b_Gt)
        
        self.calcb_Gi()
        while abs(( self.b_Ge - self.b_Gi )/ self.b_Ge ) > 0.001 : 
            self.b_Ge = min( self.b_Gi, self.b_Gt)
            self.calcb_Gi()
            self.b_Gi = min( self.b_Gi, self.b_Gt)
                        

    def calcd_Ge(self):
        """tabulka 1 (68)"""
        self.d_Ge = self.d_G2 - self.b_Ge

    def calcA_Ge(self):
        """(56)"""
        self.calcb_Ge()
        self.calcd_Ge()
        self.A_Ge = pi * self.d_Ge * self.b_Ge

    def calcQ_G0(self):
        """(57)"""
        self.calcA_Ge()
        self.calcF_G0()
        self.Q_G0 = self.F_G0 / self.A_Ge

    def calcE_G0(self):
        """(58)"""
        self.E_G0 = 1000

    def calcE_Gm(self):
        """(67)"""
        self.calcE_G0()
        self.E_Gm = 0.5 * self.E_G0

    def calcF_G0min(self):
        """(103)"""
        self.calcA_Ge()
        self.F_G0min = self.A_Ge * self.Q_A

    def calcF_Gdelta(self):
        """(105) (106)"""
        self.F_Gdelta = 0 ## neuvazujeme jiny nez I=0 zatezny stav

    def calcF_G0req(self, F_G0):
        """(107)"""
        self.F_G0 = F_G0
        self.calcF_G0min()
        self.calcF_Gdelta()
        self.F_G0req = max(self.F_G0min,self.F_Gdelta)
        return(self.F_G0req)

    def setPriruby(self, objPriruba1, objPriruba2):
        self.objPrvniPriruba = objPriruba1
        self.objDruhaPriruba = objPriruba2
