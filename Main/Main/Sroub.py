from math import *
from Soucast import *

class Sroub(Soucast):
    """description of class"""
    # vstupy
    n_B = 8
    d_B0 = 24                                                                                                       # jmenovity prumer zavitu sroubu                                [mm]
    d_Bs = 24                                                                                                       # prumer driku sroubu                                           [mm]
    d_B4 = 28
    p_t = 3                                                                                                         # stoupani zavitu                                               [mm]
    l_S = 0                                                                                                         # obr3 - delka driku sroubu                                            [mm]
    l_B = 70                                                                                                        # obr3 - delka zatizene casti sroubu                                   [mm]
    f_B0 = 300                                                                                                      # jmenovite (dovolene) napeti ve sroubu                         [MPa]
    F_B0spec = 100000
    Eps1_plus = 0
    Eps1_minus = 0
    A = 12
    mu_t = 0.2
    mu_n = 0.2
    alpha = 30
    kontrola = True
    zatizeni_sroubu = False


    # vypocty
    def calc632(self):                                                                                          # vypocet ucinne plochy prurezu sroubu
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

    def calck_B(self,d_n):
        """(B.6)"""
        #POZOR! alpha je polovicní uhel zavitu
        #(B.6):
        #self.k_B = self.p_t / (2 * pi) + self.mu_t * self.d_B0 * 0.9 / (2 * numpy.cos(self.alpha*pi/180)) + self.mu_n * d_n /2
        
        #vzorec vychazi z hodnoty v excelu 0.519*d_B0*mu_t, ve tvaru rce (B.6), bez posledního člene
        self.k_B = self.p_t / (2 * pi) + self.mu_t * self.d_B0 * 0.9 / (2 * numpy.cos(self.alpha*pi/180))

    def calcPhi_B(self):
        """(123)"""
        self.calcM_tBnom()
        self.calcF_BI()
        self.calcc_A()
        self.calcc_B()
        self.F_B = numpy.insert(self.F_BI,[0],self.F_B0max,1)
        self.F_B1 = numpy.insert(self.F_BI,[0],self.F_B0max1,1)
        self.F_BEXCEL = numpy.insert(self.F_BIEXCEL,[0],self.F_B0nom,1)
        
        self.Phi_B = (1 / (self.objSrouby.f_B0 * self.c_B)) * ((self.F_B/self.objSrouby.A_B)**2 + \
           3*(self.c_AI * self.M_tBnom*1000 / self.objSrouby.l_B)**2)**(1/2)
        self.Phi_B1 = (1 / (self.objSrouby.f_B0 * self.c_B)) * ((self.F_B1/self.objSrouby.A_B)**2 + \
           3*(self.c_AI * self.M_tBnom*1000 / (pi*self.objSrouby.d_Bs**3/16))**2)**(1/2)
        self.Phi_B2 = (1 / (self.objSrouby.f_B0 * self.c_B)) * ((self.F_B/self.objSrouby.A_B)**2 + \
           3*(self.c_AI * self.M_tBnom*1000 / (pi*self.objSrouby.d_Bs**3/16))**2)**(1/2)
        self.Phi_BEXCEL = ((1 / (self.objSrouby.f_B0 * self.c_B))) * ((self.F_BEXCEL/self.objSrouby.A_B)**2 + \
           3*(self.c_AI * self.M_tnomEXCEL / (pi*self.objSrouby.d_Bs**3/16))**2)**(1/2)

    def calcc_A(self):
        """(124)(125)(126)"""
        self.c_AI = numpy.zeros((len(self.P_I)))
        if self.objSrouby.A >= 10:
            self.c_AI[0] = 1
        elif self.objSrouby.A < 10:
            self.c_AI[0] = 4/3

    def calcc_B(self):
        """(127)"""
        self.c_B = min(1, self.objMatice.e_N * self.objMatice.f_N / 8 * self.objSrouby.d_B0 * self.objSrouby.f_B0,\
           self.objPriruba2.l_5t * self.objPriruba2.f_F / 8 * self.objSrouby.d_B0 * self.objSrouby.f_B0)
        if self.c_B < 1:
            print('Konstrukci lze zlepsit, protoze c_B < 1!')
