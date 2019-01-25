from math import *
from PrirubaSKuzelovymKrkem import *
from Soucast import *

class TesneniTyp1(Tesneni):
    """Plocha teseni, nizka tvrdost, kompozit nebo ciste kovove materialy - obrazek 3a)"""

    def calcb_Gifirst(self): 
        """(64)"""
        self.b_Gi = self.b_Gt

    def calcb_Gi(self, F_G0):
        """(65)"""
        self.calcE_Gm()
        self.calcd_Ge()
        #self.b_Gi = self.b_Gt
        self.b_Gi = ( self.e / ( pi * self.d_Ge * self.E_Gm) / ( self.geth_G0(self.objPrvniPriruba) * self.objPrvniPriruba.Z_F / self.objPrvniPriruba.E[0]
                    + self.geth_G0(self.objDruhaPriruba) * self.objDruhaPriruba.Z_F / self.objDruhaPriruba.E[0] ) + 
                    ( F_G0 / (pi * self.d_Ge * self.Q_smax))**2)**(1/2)

    def calcE_Gm(self):
        """(66)(67)"""
        self.calcE_G0()
        if self.typ == 1:
            self.E_Gm = self.E_G0
        elif self.typ == 2:
            self.E_Gm = 0.5 * self.E_G0
        else:
            print('Zkontrolujte zvoleny typ tesneni!')
            sys.exit(int(0))

    def calcd_Ge(self):
        """(68)"""
        self.d_Ge = self.d_G2 - self.b_Ge