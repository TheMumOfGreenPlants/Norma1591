from math import *
from Soucast import *

class Sroub(Soucast):
    """description of class"""

    # vstupy
    n_B = 0
    d_B0 = 0                                                                                                       # jmenovity prumer zavitu sroubu                                [mm]
    d_Bs = 0                                                                                                       # prumer driku sroubu                                           [mm]
    d_B4 = 0
    p_t = 0                                                                                                         # stoupani zavitu                                               [mm]
    l_S = 0                                                                                                         # obr3 - delka driku sroubu                                            [mm]
    l_B = 0                                                                                                        # obr3 - delka zatizene casti sroubu                                   [mm]
    f_B0 = 0                                                                                                      # jmenovite (dovolene) napeti ve sroubu                         [MPa]
    F_B0spec = 0
    Eps1_plus = 0
    Eps1_minus = 0
    A = 0
    mu_t = 0
    mu_n = 0
    alpha = 30
    kontrola = True

    # vypocty
    def VypocitejSrouby(self):
        self.calc632()
        self.calc633()

    def calc632(self):                                                                                              # vypocet ucinne plochy prurezu sroubu
        """(41)"""                                                                                                  # cislo rovnice
        self.d_Be = self.d_B0 - 0.9382 * self.p_t                                                                   # priloha A - tabA.1 !!!upraveny vzorec podle excelu!!!!
        self.A_B = ( min( self.d_Be , self.d_Bs ) )**2 * self.n_B * pi / 4                                               # ucinna plocha prurezu sroubu                                  [mm^2]

    def calc633(self):                                                                                                  # vypocet osoveho modulu pruznosti sroubu
        """(42)"""                                                                                                      # cislo rovnice
        self.l_e = self.l_B - self.l_S                                                                                  # delka zatizene casti zavitu - viz str. 13                     [mm]
        self.X_B = ( self.l_S / self.d_Bs**2 + self.l_e / self.d_Be**2 + 0.8 / self.d_B0 ) * 4 / ( self.n_B * pi )      # osovy modul pruznosti sroubu                                  [mm^-1]
  
    def calcEps(self):
        """(B.1)(B.2)"""
        self.Eps_plus = self.Eps1_plus  * (1 + 3 / (self.n_B)**(1/2) ) / 4
        self.Eps_minus = self.Eps1_minus  * (1 + 3 / (self.n_B)**(1/2) ) / 4

    def calck_B(self):
        """(B.6)"""
        #POZOR! alpha je polovicní uhel zavitu
        #(B.6):
        #self.k_B = self.p_t / (2 * pi) + self.mu_t * self.d_B0 * 0.9 / (2 * numpy.cos(self.alpha*pi/180)) + self.mu_n * d_n /2
        
        #vzorec vychazi z hodnoty v excelu 0.519*d_B0*mu_t, ve tvaru rce (B.6), bez posledního člene
        self.k_B = self.p_t / (2 * pi) + self.mu_t * self.d_B0 * 0.9 / (2 * numpy.cos(self.alpha*pi/180))

