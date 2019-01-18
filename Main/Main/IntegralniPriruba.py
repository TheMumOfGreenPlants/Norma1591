from Priruba import *

class IntegralniPriruba(Priruba):
    """description of class"""

    e_P = 31       # cast tloustky priruby radialne zatizena tlakem [mm]
    skorepina = 1   # 1 - kúzelova nebo valcova skorepina; 2 - kulova skorepina
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
        def f(skorepina):
            return {
                1 : 0.85 / cos ( self.Fi_S ),        # pro kuzelovou nebo valcovou skorepinu
                2 : 0.35 / cos ( self.Fi_S ),        # pro kulovou skorepinu
                }[skorepina]
        self.k_Q = f(skorepina)

    def calck_R(self, skorepina):
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

    def calch_Q(self):
        """(79)"""
        self.h_Q = (self.h_S)