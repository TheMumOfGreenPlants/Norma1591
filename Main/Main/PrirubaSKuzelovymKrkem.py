from Priruba import *
class PrirubaSKuzelovymKrkem(Priruba):
    """description of class"""
    e_1 = 26.5     # nejmensi tloustka steny na tenkem konci krku   [mm]
    e_2 = 28.5     # tloustka steny na silnem konci krku            [mm]
    l_H = 40       # delka krku                                     [mm]
    Fi_S = 0       # natoceni pripojne skorepiny                    [rad]
    d_1 = 91.5     # stredni prumer krku na tenci strane            [mm]
    d_2 = 93.5     # stredni prumer krku na silnejsi strane         [mm]
    e_P = 31       # cast tloustky priruby radialne zatizena tlakem [mm]
    skorepina = 1
    # !!!nutno vytvorit metodu - asi u GUI

    def sete(self):
        self.e = self.e_Ft

    def calce_E(self):
        """(17)"""
        self.calcBeta()
        self.e_E = self.e_1 * (1 + ( self.Beta - 1 ) * self.l_H / ( ( self.Beta / 3 ) * sqrt ( self.d_1 * self.e_1 ) + self.l_H ) )

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

    def calcGama(self):
        """(25)"""
        self.calce_E()
        self.calcd_E()        
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
                1 : 0.85 / cos ( self.Fi_S ),        # pro kuzelovou nabo valcovou skorepinu
                2 : 0.35 / cos ( self.Fi_S ),        # pro kulovou skorepinu
                }[skorepina]
        self.k_Q = f(skorepina)

    def calck_R(self, skorepina):
        """(33)"""
        def f(skorepina):
            return {
                1 : - 0.15 / cos ( self.Fi_S ),        # pro kuzelovou nabo valcovou skorepinu
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
        self.h_Q = (self.h_S * self.k_Q + self.h_T * (2 * self.d_F * self.e_P / self.d_E**2 - 0.5 * tan(self.Fi_S))) * (self.d_E / d_Ge)**2
        self.h_G = (self.d_3e - d_Ge) / 2
        self.h_H = (self.d_3e - self.d_E) / 2
        self.h_L = 0