from math import *

class Sroub(object):
    """description of class"""
    # vstupy
    d_B0 = 24                                                                                                       # jmenovity prumer zavitu sroubu                                [mm]
    d_Bs = 24                                                                                                       # prumer driku sroubu                                           [mm]
    p_t = 3                                                                                                         # stoupani zavitu                                               [mm]
    l_S = 0                                                                                                         # delka driku sroubu                                            [mm]
    l_B = 70                                                                                                        # delka zatizene casti sroubu                                   [mm]

    # vypocty
    def calcA_B(self,n_B):                                                                                          # vypocet ucinne plochy prurezu sroubu
        """(41)"""                                                                                                  # cislo rovnice
        self.d_Be = self.d_B0 - 0.9382 * self.p_t                                                                   # priloha A - tabA.1 !!!upraveny vzorec podle excelu!!!!
        self.A_B = ( min( self.d_Be , self.d_Bs ) )**2 * n_B * pi / 4                                               # ucinna plocha prurezu sroubu                                  [mm^2]

    def calcX_B(self,n_B):                                                                                          # vypocet osoveho modulu pruznosti sroubu
        """(42)"""                                                                                                  # cislo rovnice
        self.l_e = self.l_B - self.l_S                                                                              # delka zatizene casti zavitu - viz str. 13                     [mm]
        self.X_B = ( self.l_S / self.d_BS**2 + self.l_e / self.d_Be**2 + 0.8 / self.d_B0 ) * 4 / ( n_B * pi )       # osovy modul pruznosti sroubu                                  [mm^-1]

    def calcF_R0(self):                                                                                             # vypocet pridavnych vnejsich zatizeni
        """(96)"""                                                                                                  # cislo rovnice
        self.F_R0 = 0                                                                                               # bez pridavnych vnejsich zatizeni                              [N] 

    def calcF_B0req(self,F_G0req):                                                                                  # vypocet pozadovane sily ve sroubech - MONTAZNI STAV
        """(108)"""                                                                                                 # cislo rovnice
        self.calcF_R0()                                                                                             # k vypoctu potrebujeme znat F_R0
        self.F_B0req = F_G0req + self.F_R0                                                                          # pozadovana sila ve sroubech                                   [N]

    def calcF_B0nom(self, F_G0req):                                                                                 # vypocet skutecne sily pusobici na srouby - MONTAZNI STAV
        """(112)"""                                                                                                 # cislo rovnice
        self.calcF_B0req(F_G0req)                                                                                   # k vypoctu potrebujeme znat F_B0req
        self.F_B0nom = self.F_B0req / (1 - 0) # 0 = epsilon- zjenodušno!                                            # skutecna sila pusobici na srouby                              [N]
    
    def calcPreload(self, F_G0req):                                                                                 # vypocet predpeti ve sroubu
        self.calcF_B0nom(F_G0req)                                                                                   # k vypoctu potrebujeme znat F_B0nom
        self.Preload = self.F_B0nom / self.A_B                                                                      # predpeti ve sroubu                                            [MPa]