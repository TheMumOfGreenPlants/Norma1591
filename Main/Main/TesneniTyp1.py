from math import *
from Soucast import *
from Tesneni import *

class TesneniTyp1(Tesneni):
    """Plocha teseni, nizka tvrdost, kompozit nebo ciste kovove materialy - obrazek 3a)"""

    druh = 0  # 1 - ploche kovove kruhove tesneni s pravouhlym prurezem; 2 - pro nekovova plocha tesneni

    def calcb_Gifirst(self): 
        """(64)"""
        self.b_Gi = self.b_Gt

    def calcb_Gi(self,obj1,obj2,F_G0):
        """(65)"""
        self.b_Gi = ( self.e / ( pi * self.d_Ge * self.E_Gm)
                    / ( obj1.h_G0 * obj1.Z_F / obj1.E[0] + obj2.h_G0 * obj2.Z_F / obj2.E[0] ) + 
                    ( F_G0 / (pi * self.d_Ge * self.Q_smax))**2)**(1/2)



    def calcE_Gm(self,F_G0):
        """(66)(67)"""
        if self.druh == 1:
            self.E_Gm = self.E_G0
        elif self.druh == 2:
            self.E_Gm = 0.5 * self.E_G0
        else:
            print('Zkontrolujte zvoleny druh tesneni!')
            sys.exit(int(0))

    def calcd_Ge(self):
        """(68)"""
        self.d_Ge = self.d_G2 - self.b_Ge