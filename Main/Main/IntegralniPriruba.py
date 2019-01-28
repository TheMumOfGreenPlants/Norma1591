from Priruba import *
from numpy import *

class IntegralniPriruba(Priruba):
    """description of class"""

    e_P = 31        # cast tloustky priruby radialne zatizena tlakem [mm]
    skorepina = 1   # 1 - kuzelova nebo valcova skorepina; 2 - kulova skorepina
    # !!!nutno vytvorit metodu - asi u GUI

    def calcbde_FL(self):
        """(7)(8)(9)(10)"""
        self.b_F = (self.d_4 - self.d_0) / 2 - self.d_5e
        self.b_L = 0
        self.d_L = 0
        self.e_L = 0
        self.d_F = (self.d_4 + self.d_0)/2
        #self.e_F = 2 * self.A_F / (self.d_4 - self.d_0)

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
        self.h_H = (self.d_3e - self.d_Ed_E) / 2
        self.h_L = 0

    def calcPhi_F(self):
        """(129)"""
        p.calce_D()
        self.calcW_F()
        self.Phi_F = abs(self.F_G * self.h_G + self.F_QI * (self.h_H - self.h_P) + self.F_RI * self.h_H) / self.W_F

    def calcW_F(self):
        """(130)"""
        self.W_F = (pi/4)*(self.f_F * 2 * self.b_F * self.e**2 * (1 + 2 * self.Psi_opt * self.Psi_Z))

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
                1 : (1.33 * (1 - 0.75 * (0.5 * self.delta_Q + self.delta_R)**2 * (1 - (0.75 * self.delta_Q**2 + 1 * self.delta_R**2))))**(1/2),        # pro kuzelovou nebo valcovou skorepinu
                2 : (1.33 * (1 - 0.75 * (0.5 * self.delta_Q + self.delta_R)**2 * (1 - (0.25 * self.delta_Q**2 + 3 * self.delta_R**2))))**(1/2),        # pro kulovou skorepinu
                }[skorepina]
        self.c_M = f(self.skorepina)

    def calcc_S(self):
        """(135)"""
        def f(skorepina):
            j_s = numpy.asarray([-1,1])
            return {
                1 : (pi /4) * (1 - 0.75 * (0.5 * self.delta_Q + self.delta_R)**2)**(1/2) + j_s * (0.5 * self.delta_R - 0.75 * self.delta_Q),        # pro kuzelovou nebo valcovou skorepinu
                2 : (pi /4) * (1 - 0.75 * (0.5 * self.delta_Q + self.delta_R)**2)**(1/2) + j_s * (1.5 * self.delta_R - 0.25 * self.delta_Q),        # pro kulovou skorepinu
                }[skorepina]
        self.c_Sminus = f(self.skorepina)


    def calcj_M(self):
        """(136)"""
        self.j_M = 2

