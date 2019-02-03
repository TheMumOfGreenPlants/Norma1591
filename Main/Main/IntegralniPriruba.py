from Priruba import *
from numpy import *
import sys

class IntegralniPriruba(Priruba):
    """description of class"""
    def __init__(self,skorepina_vstup,krk_vstup):
        self.skorepina = skorepina_vstup
        self.krk = krk_vstup

    e_P = 31        # cast tloustky priruby radialne zatizena tlakem [mm]
    skorepina = 1   # 1 - kuzelova nebo valcova skorepina; 2 - kulova skorepina
    # !!!nutno vytvorit metodu - asi u GUI
    e_1 = 26.5      # nejmensi tloustka steny na tenkem konci krku   [mm]
    e_2 = 28.5      # tloustka steny na silnem konci krku            [mm]
    l_H = 40        # delka krku                                     [mm]
    Fi_S = 0        # natoceni pripojne skorepiny                    [rad]
    d_1 = 91.5      # stredni prumer krku na tenci strane            [mm]
    d_2 = 93.5      # stredni prumer krku na silnejsi strane         [mm]
    krk = 1         # 1 - ano, 2 - ne
    j_S = numpy.asarray([-1,1])

    def calcbde_FL(self):
        """(7)(8)(9)(10)"""
        self.b_F = (self.d_4 - self.d_0) / 2 - self.d_5e
        self.b_L = 0
        self.d_L = 0
        self.e_L = 0
        self.d_F = (self.d_4 + self.d_0)/2
        #self.e_F = 2 * self.A_F / (self.d_4 - self.d_0)




    def calce_E(self):
        """(17)(21)"""
        self.calcBeta()
        def switche_E(self):
            return {
                1 : self.e_1 * (1 + ( self.Beta - 1 ) * self.l_H / ( ( self.Beta / 3 ) * sqrt ( self.d_1 * self.e_1 ) + self.l_H ) ),        # priruba s krkem
                2 : self.e_S,        # Priruba bez krku. Vzorec neplati v pripade, kdy je hrdlo pripojeno ke stredovemu otvoru zaslepovaci priruby. Pro tento pripad plati (23)!"""
                }[self.krk]
        self.e_E = switche_E()

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
        # (22)
        # Vzorec neplati v pripade, kdy je hrdlo pripojeno ke stredovemu otvoru zaslepovaci priruby. Pro tento pripad plati (24)!
        
        self.d_E = self.d_S


    def calcGama(self):
        """(25)"""
        self.calce_E()
        self.calcd_E()
        self.calcd_5e()
        self.calcbde_FL()
        self.Gama = self.e_E * self.d_F / ( self.b_F * self.d_E * cos( self.Fi_S ) )

    def calcTheta(self):
        """(26)"""
        self.Theta = 0.55 * cos( self.Fi_S ) * sqrt( self.d_E * self.e_E ) / self.e_F * 1      # zkontrolovat vzorec!

    def calcLambda(self):
        """(27)"""
        self.Lambda = 1 - self.e_P / self.e_F

    def calcc_F(self):
        """(28)"""
        self.calcGama()
        self.calcTheta()
        self.calcLambda()
        self.c_F = ( 1 + self.Gama * self.Theta ) / ( 1 + self.Gama * self.Theta * ( 4 *( 1 - 3 * self.Lambda + 3 * self.Lambda**2) + \
            6 * ( 1 - 2 * self.Lambda) * self.Theta + 6 * self.Theta**2 ) + 3 * self.Gama**2 * self.Theta**4)

    def calch_S(self):
        """(29)"""
        self.h_S = 1.1 * self.e_F * sqrt( self.e_E / self.d_E ) * ( 1 - 2 * self.Lambda + self.Theta ) / ( 1 + self.Gama * self.Theta )

    def calch_T(self):
        """(30)"""
        self.h_T = self.e_F * (1 - 2 * self.Lambda - self.Gama * self.Theta**2) / ( 1 + self.Gama * self.Theta )

    def calch_R(self):
        """(31)"""
        self.calcGama()
        self.calcTheta()
        self.calcLambda()
        self.calch_S()
        self.calch_T()
        self.calck_R(self.skorepina)
        self.h_R = self.h_S * self.k_R - self.h_T * 0.5 * tan( self.Fi_S )

    def calck_Q(self, skorepina):
        """(32)"""
        def f(self):
            return {
                1 : 0.85 / cos ( self.Fi_S ),        # pro kuzelovou nebo valcovou skorepinu
                2 : 0.35 / cos ( self.Fi_S ),        # pro kulovou skorepinu
                }[self.skorepina]
        self.k_Q = f(skorepina)

    def calck_R(self):
        """(33)"""
        def f(skorepina):
            return {
                1 : - 0.15 / cos ( self.Fi_S ),        # pro kuzelovou nebo valcovou skorepinu
                2 : - 0.65 / cos ( self.Fi_S ),        # pro kulovou skorepinu
                }[skorepina]
        self.k_R = f(skorepina)

    def calcZ_FL(self):
        """(34) (35)"""
        self.calcc_F()
        self.Z_F = 3 * self.d_F * self.c_F / ( pi * self.b_F * self.e_F**3 )
        self.Z_L = 0


    def calch_QGHL(self, d_Ge):
        """(79)(81)(82)(83)"""
        self.calck_Q(self.skorepina)
        self.h_Q = (self.h_S * self.k_Q + self.h_T * (2 * self.d_F * self.e_P / self.d_E**2 - 0.5 * tan(self.Fi_S))) * (self.d_E / objTesneni.d_Ge)**2
        self.h_G = (self.d_3e - d_Ge) / 2
        self.h_G = (self.d_3e - objTesneni.d_Ge) / 2
        self.h_H = (self.d_3e - self.d_E) / 2
        self.h_L = 0

    def calcf_E(self):
        """(131)"""
        self.f_E = min(self.f_F,self.f_S)

    def calcdelta_Q(self):
        """(132)"""
        self.delta_Q = self.objZatizeni.P_I * self.d_E / (self.f_E * 2 * self.e_D * cos(self.Fi_S))

    def calcdelta_R(self):
        """(133)"""
        self.delta_R = self.objZatizeni.F_RI / (self.f_E * pi * self.d_E * self.e_D * cos(self.Fi_S))

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
        def f(skorepina):
            return {
                1 : (pi /4) * (1 - 0.75 * (0.5 * self.delta_Q + self.delta_R)**2)**(1/2) + self.j_S * (0.5 * self.delta_R - 0.75 * self.delta_Q),        # pro kuzelovou nebo valcovou skorepinu
                2 : (pi /4) * (1 - 0.75 * (0.5 * self.delta_Q + self.delta_R)**2)**(1/2) + self.j_S * (1.5 * self.delta_R - 0.25 * self.delta_Q),        # pro kulovou skorepinu
                }[skorepina]
        self.c_S = f(self.skorepina)


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
        if (self.Psi_max < -1) or (self.Psi_min > 1):
            print("List priruby je pretizen!")
            sys.exit(int(0))

    def calcPsi_Z(self):
        """tabulka 2 - urcovani Psi_Z"""
        if self.j_M == 1:
            if self.Psi_max <= self.Psi_opt:
                self.k_M = 1
                self.Psi_Z = self.Psi_max
            elif (self.Psi_0 <= self.Psi_opt) and (self.Psi_opt < self.Psi_max):
                self.k_M = 1
                self.Psi_Z = self.Psi_opt
            elif self.Psi_opt < self.Psi_0:
                self.k_M = 1 # POZOR!-nedoděláno norma: k_M < +1
                self.Psi_Z = getPsi(-1,self.k_M,1)
        elif self.j_M == -1:
            if self.Psi_opt <= self.Psi_min:
                self.k_M = -1
                self.Psi_Z = self.Psi_min
            elif (self.Psi_min < self.Psi_opt) and (self.Psi_opt <= self.Psi_0):
                self.k_M = -1
                self.Psi_Z = self.Psi_opt
            elif self.Psi_0 < self.Psi_opt:
                self.k_M = -1 # POZOR!-nedoděláno norma: k_M > -1
                self.Psi_Z = getPsi(1,self.k_M,1)

    def calcW_F(self):
        """(130)"""
        self.calcPsi_Z
        self.W_F = (pi / 4) * (self.f_F * 2 * self.b_F * self.e_F**2 * (1 + 2 *self.Psi_opt *self.Psi_Z - self.Psi_Z**2) + self.f_E * self.d_E * self.e_D**2 * self.c_M * self.j_M * self.k_M)

    def calcPhi_F(self):
        """(129)"""
        self.Phi_F = abs(F_G * self.h_G + F_Q * (self.h_H - self.h_P) + F_R * self.h_H) / self.W_F