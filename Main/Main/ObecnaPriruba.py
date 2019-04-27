from Priruba import *
from numpy import *
import sys

class ObecnaPriruba(Priruba):
    """description of class"""

    def __init__(self, krk_volba, skorepina_volba, diry_volba):
        self.krk = krk_volba
        self.skorepina = skorepina_volba
        self.diry = diry_volba


# VSTUPNI PARAMETRY
    e_P = 0        # cast tloustky priruby radialne zatizena tlakem [mm]
    e_1 = 0      # nejmensi tloustka steny na tenkem konci krku   [mm]
    e_2 = 0      # tloustka steny na silnem konci krku            [mm]
    l_H = 0        # delka krku                                     [mm]
    d_1 = 0      # stredni prumer krku na tenci strane            [mm]
    d_2 = 0      # stredni prumer krku na silnejsi strane         [mm]
# KONEC - VSTUPNI PARAMETRY

    #j_S = numpy.asarray([-1,1])

    def beforecalc(self):
        # vsechny typy krome 4 a 5:
        if not((self.krk == 0) & (self.skorepina == 2) & (self.diry == 2)) | ((self.krk == 1) & (self.skorepina == 3) & (self.diry == 1)):
            self.e_P = self.e_F
        # typ 2:
        if ((self.krk == 1) & (self.skorepina == 1) & (self.diry == 2)):
            self.e_S = self.e_1
        # typ 7:
        if ((self.krk == 0) & (self.skorepina == 1) & (self.diry == 1)):
            self.e_S = self.e_1
            self.e_2 = self.e_1
        # typ 11:
        if ((self.krk == 1) & (self.skorepina == 1) & (self.diry == 1)):
            self.e_S = self.e_1

    def calc622(self):
        """(5) (7)(8)(9) *(10)*"""
        if self.diry == 2:
            self.d_5 = self.d_5t * self.l_5t / self.e_Fb
        self.calc6221()
        self.b_F = (self.d_4 - self.d_0) / 2 - self.d_5e
        self.b_L = 0
        self.d_L = 0
        self.e_L = 0
        self.d_F = (self.d_4 + self.d_0)/2
        #self.e_F = 2 * self.A_F / (self.d_4 - self.d_0)

    def calc623(self):
        if self.krk == 1:       # pro priruby S krkem
            self.calc6231()
        elif self.krk == 0:     # pro priruby BEZ krku
            self.calc6232()

    def calc6231(self):
        """(19)(17)(18)(20)"""
        self.Beta = self.e_2 / self.e_1
        self.e_E = self.e_1 * (1 + (self.Beta - 1) * self.l_H / ((self.Beta / 3) * sqrt (self.d_1 * self.e_1) + self.l_H))
        self.e_D = self.e_1 * ( 1 + ( self.Beta - 1 ) * self.l_H / ( (self.Beta / 3)**4 * (self.d_1 * self.e_1)**2 + self.l_H**4 )**(1/4))
        self.d_E = (min((self.d_1 - self.e_1 + self.e_E), (self.d_2 + self.e_2 - self.e_E)) + max((self.d_1 + self.e_1 - self.e_E), (self.d_2 - self.e_2 + self.e_E))) / 2

    def calc6232(self):
        """(21)(22)"""
        # Vzorec neplati v pripade, kdy je hrdlo pripojeno ke stredovemu otvoru zaslepovaci priruby. Pro tento pripad plati (23)!
        self.e_E = self.e_S
        # Vzorec neplati v pripade, kdy je hrdlo pripojeno ke stredovemu otvoru zaslepovaci priruby. Pro tento pripad plati (24)!
        self.d_E = self.d_S

    def calc624(self):
        """(25)(26)(27)(28)(29)(30)(31)(32)(33)(34)(35)"""
        self.Gama = self.e_E * self.d_F / ( self.b_F * self.d_E * cos( self.Fi_S ) )
        self.Theta = 0.55 * cos( self.Fi_S ) * sqrt( self.d_E * self.e_E ) / self.e_F * 1      # zkontrolovat vzorec!
        self.Lambda = 1 - self.e_P / self.e_F
        self.c_F = ( 1 + self.Gama * self.Theta ) / ( 1 + self.Gama * self.Theta * ( 4 *( 1 - 3 * self.Lambda + 3 * self.Lambda**2) + \
            6 * ( 1 - 2 * self.Lambda) * self.Theta + 6 * self.Theta**2 ) + 3 * self.Gama**2 * self.Theta**4)
        self.h_S = 1.1 * self.e_F * sqrt( self.e_E / self.d_E ) * ( 1 - 2 * self.Lambda + self.Theta ) / ( 1 + self.Gama * self.Theta )
        self.h_T = self.e_F * (1 - 2 * self.Lambda - self.Gama * self.Theta**2) / ( 1 + self.Gama * self.Theta )
        def calck_Q(self):
            return {
                1 : 0.85 / cos ( self.Fi_S ),        # pro valcovou skorepinu
                2 : 0.85 / cos ( self.Fi_S ),        # pro kuzelovou skorepinu
                3 : 0.35 / cos ( self.Fi_S ),        # pro kulovou skorepinu
                }[self.skorepina]
        def calck_R(self):
            return {
                1 : - 0.15 / cos ( self.Fi_S ),         # pro valcovou skorepinu
                2 : - 0.15 / cos ( self.Fi_S ),         # pro kuzelovou skorepinu
                3 : - 0.65 / cos ( self.Fi_S ),         # pro kulovou skorepinu
                }[self.skorepina]
        self.k_Q = calck_Q(self)
        self.k_R = calck_R(self)
        self.h_R = self.h_S * self.k_R - self.h_T * 0.5 * tan( self.Fi_S )
        self.Z_F = 3 * self.d_F * self.c_F / ( pi * self.b_F * self.e_F**3 )
        self.Z_L = 0

    def calc645(self,d_Ge):
        self.calch_P(d_Ge)
        #self.calck_Q(self.skorepina)
        self.calch_Q(d_Ge)
        self.calch_GHL(d_Ge)

    def calch_GHL(self,d_Ge):
        """(81)(82)(83)"""
        self.h_G = (self.d_3e - d_Ge) / 2
        self.h_H = (self.d_3e - self.d_E) / 2
        self.h_L = 0

    def calch_Q(self,d_Ge):
        """(79)"""
        self.h_Q = (self.h_S * self.k_Q + self.h_T * (2 * self.d_F * self.e_P / self.d_E**2 - 0.5 * tan(self.Fi_S))) * (self.d_E / d_Ge)**2



    def calch_G0(self, d_Ge):
        """(59)"""
        self.h_G0 = (self.d_3e - d_Ge) / 2

    def calcf_E(self):
        """(131)"""
        if self.l_H == 0:
            self.e_D = self.e_1
        self.f_E = min(self.f_F,self.f_S)

    def calcdelta_QR(self,P_I,F_RI):
        """(132)"""
        self.delta_Q = P_I * self.d_E / (self.f_E * 2 * self.e_D * cos(self.Fi_S))
        self.delta_R = F_RI / (self.f_E * pi * self.d_E * self.e_D * cos(self.Fi_S))

    def calcc_M(self):
        """(134)"""
        def f(skorepina):
            return {
                1 : (1.33 * (1 - 0.75 * (0.5 * self.delta_Q + self.delta_R)**2) * (1 - (0.75 * self.delta_Q**2 + 1 * self.delta_R**2)))**(1/2),        # pro kuzelovou nebo valcovou skorepinu
                2 : (1.33 * (1 - 0.75 * (0.5 * self.delta_Q + self.delta_R)**2) * (1 - (0.25 * self.delta_Q**2 + 3 * self.delta_R**2)))**(1/2),        # pro kulovou skorepinu
                }[skorepina]
        self.c_M = f(self.skorepina)
        if isinstance(self.c_M,complex):
            print("Krk je pretizeni!")
            sys.exit(int(0))

    def calcc_S(self):
        """(135)"""
        def f(skorepina,j_S):
            return {
                1 : (pi /4) * (1 - 0.75 * (0.5 * self.delta_Q + self.delta_R)**2)**(1/2) + j_S * (0.5 * self.delta_R - 0.75 * self.delta_Q),        # pro kuzelovou nebo valcovou skorepinu
                2 : (pi /4) * (1 - 0.75 * (0.5 * self.delta_Q + self.delta_R)**2)**(1/2) + j_S * (1.5 * self.delta_R - 0.25 * self.delta_Q),        # pro kulovou skorepinu
                }[skorepina]
        self.c_S = [f(self.skorepina,1),f(self.skorepina,-1)]


    def calcj_M(self,F_G,F_Q,F_R):
        """(136)"""
        self.j_M = (F_G * self.h_G + F_Q * (self.h_H - self.h_P) + F_R * self.h_H) / abs(F_G * self.h_G + F_Q * (self.h_H - self.h_P) + F_R * self.h_H)

    def getPsi(self,j_S,k_M,k_S):
        """(140)"""
        Psi = (self.f_E * self.d_E * self.e_D * cos(self.Fi_S)) / (self.f_F * 2 * self.b_F * self.e_F) \
            * ((0.5 * self.delta_Q + self.delta_R) * tan(self.Fi_S) - self.delta_Q * 2 * self.e_P / self.d_E \
            + j_S * k_S * ((self.e_D * self.c_M * self.c_S * (1 + j_S * k_M)) / (self.d_E * (cos(self.Fi_S))**3))**(1/2))
        return Psi

    def calcPsi(self):
        """(141)(142)(143)(144)"""
        self.Psi_opt = self.j_M * (2 * self.e_P / self.e_F - 1)
        self.Psi_0 = self.getPsi(0,0,0)
        self.Psi_max = self.getPsi(1,1,1)
        self.Psi_min = self.getPsi(-1,-1,1)
        if ((self.Psi_max < -1).any or (self.Psi_min > 1).any) == True:
            print("List priruby je pretizen!")
            sys.exit(int(0))

    def calcPsi_Z(self):
        """tabulka 2 - urcovani Psi_Z"""
        k_M = zeros([2,len(self.j_M)])
        Psi_Z = zeros([2,len(self.j_M)])
        for n in range(0,2):
            for i,j_M in enumerate(self.j_M):
                if j_M == 1:
                    if self.Psi_max[n,i] <= self.Psi_opt[i]:
                        k_M[n,i] = 1
                        Psi_Z[n,i] = self.Psi_max[n,i]
                    elif (self.Psi_0[n,i] <= self.Psi_opt[i]) and (self.Psi_opt[i] < self.Psi_max[n,i]):
                        k_M[n,i] = 1
                        Psi_Z[n,i] = self.Psi_opt[i]
                    elif self.Psi_opt[i] < self.Psi_0[n,i]:
                        k_M[n,i] = 1 # POZOR!-nedoděláno norma: k_M < +1
                        Psi_Z[n,i] = getPsi(-1,k_M[n,i],1)
                elif j_M == -1:
                    if self.Psi_opt[i] <= self.Psi_min[n,i]:
                        k_M[n,i] = -1
                        Psi_Z[n,i] = self.Psi_min[n,i]
                    elif (self.Psi_min[n,i] < self.Psi_opt[i]) and (self.Psi_opt[i] <= self.Psi_0[n,i]):
                        k_M[n,i] = -1
                        Psi_Z[n,i] = self.Psi_opt[i]
                    elif self.Psi_0[n,i] < self.Psi_opt[i]:
                        k_M[n,i] = -1 # POZOR!-nedoděláno norma: k_M > -1
                        Psi_Z[n,i] = getPsi(1,k_M[n,i],1)
        self.k_M = zeros(len(self.j_M))
        self.Psi_Z = zeros(len(self.j_M))
        for i,val in enumerate(k_M[0]):
            if val >= k_M[1,i]:
                self.k_M[i] = val
                self.Psi_Z[i] = Psi_Z[0,i]
            else:
                self.k_M[i] = k_M[1,i]
                self.Psi_Z[i] = Psi_Z[1,i]

    def calcW_F(self):
        """(130)"""
        self.W_F = (pi / 4) * (self.f_F * 2 * self.b_F * self.e_F**2 * (1 + 2 *self.Psi_opt *self.Psi_Z - self.Psi_Z**2) + self.f_E * self.d_E * self.e_D**2 * self.c_M * self.j_M * self.k_M)

    def calc8456(self,P,F_G,F_Q,F_R,F_B,Tesneni):
        """(129)"""
        self.calcf_E()
        self.calcdelta_QR(P,F_R)
        self.calcc_M()
        self.calcc_S()
        self.calcj_M(F_G,F_Q,F_R)
        self.calcPsi()
        self.calcPsi_Z()
        self.calcW_F()
        self.Phi_F = abs(F_G * self.h_G + F_Q * (self.h_H - self.h_P) + F_R * self.h_H) / self.W_F
        return self.Phi_F
    
    def calc4211(self):
        return 0.2 <= (self.b_F / self.e_F)
        
    def calc4212(self):
        return (self.b_F / self.e_F) <= 5

    def calc421(self):
        podm11 = self.calc4211()
        podm12 = self.calc4212()
        return podm11 & podm12