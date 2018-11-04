from math import pi, acos
from Priruba import *
from Sroub import *
from Tesneni import *
from Soucast import *
from Matice import *

class Zatizeni(object):
    """description of class"""
    P_I = 0
    F_XI = 0
    F_YI = 0
    F_ZI = 0
    M_XI = 0
    M_YI = 0
    M_ZI = 0

    T_0 = 0
    T_BI = 0
    T_1FI = 0
    T_2FI = 0
    T_1LI = 0
    T_2LI = 0
    T_1WI = 0
    T_2WI = 0
    T_GI = 0

    alfa_BI = 0
    alfa_1FI = 0
    alfa_2FI = 0
    alfa_1LI = 0
    alfa_2LI = 0
    alfa_1WI = 0
    alfa_2WI = 0
    alfa_GI = 0

    F_GdeltaI = []
    F_GI = []
    F_BI = []
    c_AI = []
    F_R0 = None
    Y_R0 = None
    Y_G0 = None
    N_R = 1
        
    def setall(self, priruba1, priruba2, srouby, tesneni, matice):
        self.objPriruba1 = priruba1
        self.objPriruba2 = priruba2
        self.objSrouby = srouby
        self.objTesneni = tesneni
        self.objMatice = matice

    def calcA_Q(self):
        """(90)"""
        self.A_Q = (pi * self.objTesneni.d_Ge**2) /4     #d_Ge==d_Gi????

    def calcF_QI(self,i):
        """(91)"""
        if i == 0:
            self.calcA_Q()
            self.F_QI = []
        self.F_QI.append(self.A_Q * self.P_I[i])

    def calcFM(self, varianta):
        """(92)(93)(94)(95)(96)"""
        self.F_AI = F_ZI
        self.F_LI = (F_XI**2 + F_YI**2)**(1/2)
        self.M_AI = (M_XI**2 + M_YI**2)**(1/2)
        self.M_TGI = M_ZI
        def f(var):
            return {
                a : F_AI + (4 / self.objPriruba1.d_3e) * M_AI,
                b : F_AI - (4 / self.objPriruba1.d_3e) * M_AI,
                }[var]
        self.F_RI = f(varianta)
        if self.F_R0 == None:
            self.F_R0 = self.F_RI

    def calcdeltaU_TI(self):
        """(97)""" # !! chybi podlozky
        def part(soucast):
            return soucast.e * soucast.alfa * (soucast.T - self.T_0)
        self.deltaU_TI = self.objSrouby.l_B * self.objSrouby.alfa * (self.objSrouby.T - self.T_0) \
            - part(self.objPriruba1 ) - part(self.objPriruba2) - part(self.objTesneni)

    def conditionl_B(self, e_Ft1, e_Ft2, e_L1, e_L2, e_W1, e_W2, e_Gt, l_B):         
        """(98)"""
        sum = e_Ft1 + e_Ft2 + e_L1 + e_L2 + e_W1 + e_W2 + e_Gt                         
        return sum == l_B

    def calcY(self):
        """(99)(100)(101)(102)"""
        self.objTesneni.calcX_G()
        self.objPriruba1.calch_P(self.objTesneni)
        self.objPriruba2.calch_P(self.objTesneni)
        self.objPriruba1.calch_R()
        self.objPriruba2.calch_R()
        self.Y_BI = self.objPriruba1.Z_L * self.objPriruba1.h_L**2 / self.objPriruba1.E \
            + self.objPriruba2.Z_L * self.objPriruba2.h_L**2 / self.objPriruba2.E \
            + self.objSrouby.X_B / self.objSrouby.E           ##  BEZ PODLOZEK!!!!
        self.Y_GI = self.objPriruba1.Z_F * self.objPriruba1.h_G**2 / self.objPriruba1.E \
            + self.objPriruba2.Z_F * self.objPriruba2.h_G**2 / self.objPriruba2.E \
            + self.Y_BI + self.objTesneni.X_G/ self.objTesneni.E
        self.Y_QI = self.objPriruba1.Z_F * self.objPriruba1.h_G * (self.objPriruba1.h_H - self.objPriruba1.h_P + self.objPriruba1.h_Q)/self.objPriruba1.E \
            + self.objPriruba2.Z_F * self.objPriruba2.h_G * (self.objPriruba2.h_H - self.objPriruba2.h_P + self.objPriruba2.h_Q)/self.objPriruba2.E \
            + self.Y_BI
        self.Y_RI = self.objPriruba1.Z_F * self.objPriruba1.h_G * (self.objPriruba1.h_H + self.objPriruba1.h_R)/self.objPriruba1.E \
            + self.objPriruba2.Z_F * self.objPriruba2.h_G * (self.objPriruba2.h_H + self.objPriruba2.h_R)/self.objPriruba2.E \
            + self.Y_BI
        self.isfirst = False
        if self.Y_R0 == None:
            self.Y_R0 = self.Y_RI
            self.isfrst = True
        if self.Y_G0 == None:
            self.Y_G0 = self.Y_GI
            self.isfrst = True


    def calcF_G0min(self):
        """(103)"""
        self.F_G0min = self.objTesneni.A_Ge * self.objTesneni.Q_A

    def calcF_GImin(self):
        """(104)"""
        self.calcF_QI()
        self.calcFM()
        self.calcd_Gt()
        self.F_GImin = max(self.objTesneni.A_Ge * self.objTesneni.Q_sminLI,-(self.F_QI + self.F_RI),\
                        self.F_LI/self.objTesneni.mu_G + (2 * self.M_TGI) / (self.objTesneni.mu_G * self.objTesneni.d_Gt) - (2 * self.M_AI) / self.objTesneni.d_Gt)

    def calcdeltae_Gc(self):
        """(F.3)"""
        self.calcP_QR()
        self.deltae_Gc = self.objTesneni.K * self.Y_GI * self.objTesneni.deltae_Gc_test

    def calcF_GdeltaI(self):
        """(106) 1/2"""
        if self.isfirst == False:
            self.F_GdeltaI.append((self.F_GImin * self.Y_GI + self.F_QI * self.Y_QI + (self.F_RI * self.Y_RI - self.F_R0 * self.Y_R0) \
                + self.deltaU_TI + self.deltae_Gc + (self.objTesneni.e_G - self.objTesneni.e_GA))/ self.Y_G0)

    def calcF_Gdelta(self):
        """(106) 2/2"""
        for i in self.P_I:
            self.calcF_QI(i)
        self.F_Gdelta = max(self.F_GdeltaI)

    def calcF_G0req(self):
        """(107)"""
        self.calcF_Gdelta()
        self.F_G0req = max(self.F_G0min,self.F_Gdelta)

    def calcF_B0req(self):
        """(108)"""
        self.calcF_G0req()
        self.calcF_B0req = self.F_G0req + self.F_R0

    def calcF_B0minmax(self):
        """(112)(123)"""
        self.calcF_B0av()
        self.FB_0min = self.F_B0av * (1 - self.objSrouby.Eps_minus)
        self.F_B0max = self.F_B0av * (1 - self.objSrouby.Eps_plus)

    def calcF_B0av(self):
        """(B.3)"""
        self.F_B0av = min(self.objSrouby.A_B * self.objSrouby.f_B0, self.objSrouby.n_B *200000 )

    def calcF_B0max(self):
        """(117)"""
        self.F_B0max = self.objSrouby.F_B0nom * ( 1 + self.objSrouby.Eps_plus )

    def calcF_G0max(self):
        """(118)"""
        self.F_G0max = self.F_B0max - self.F_R0

    def calcF_G0d(self):
        """(119)"""
        self.F_G0d = max(self.F_Gdelta , (2/3) * (1 - 10 / self.N_R) * self.F_B0max - self.F_R0)

    def calcF_GI(self):
        """(121)"""
        self.calcF_Gd()
        self.F_GI.append((self.F_G0d * self.Y_G0 - (self.F_QI * self.Y_QI + (self.F_RI * self.Y_RI - self.F_R0 * self.Y_R0) \
                + self.deltaU_TI ) - self.deltae_Gc + (self.objTesneni.e_G - self.objTesneni.e_GA)) / self.Y_GI)

    def calcF_BI(self):
        """(122)"""
        self.calcF_GI()
        self.F_BI.append(self.F_GI + (self.F_QI + self.F_RI))

    def calcc_A(self):
        """(124)(125)(126)"""
        if len(self.c_AI) == 0 and self.A >= 10:
            self.c_AI.append(1)
        elif len(self.c_AI) == 0 and self.A < 10:
            self.c_AI.append(4/3)
        else:
            self.c_AI.append(0)

    def calcc_B(self):
        """(127)"""
        self.c_B = min(1, self.objMatice.e_N * self.objMatice.f_N / 8 * self.objSrouby.d_B0 * self.objSrouby.f_B0,\
           self.objPriruba2.l_5t * self.objPriruba2.f_F / 8 * self.objSrouby.d_B0 * self.objSrouby.f_B0)

    def calcPhi_B(self):
        """(123)"""
        self.calcc_A()
        self.calcc_B()
        self.calcM_tBI()
        self.Phi_B = (1 / (self.F_BI * self.objSrouby.c_B)) * ((self.F_BI/self.objSrouby.A_B)**2 + \
           3*(self.c_AI * self.M_tBI / self.objSrouby.l_B)**2)**(1/2)

    def calcM_tBI(self):
        """(B.4)"""
        self.calck_B()
        self.M_tBI = self.k_B * self.objSrouby.F_B0nom / self.objSrouby.n_B

    def calck_B(self):
        """(B.6)(B.7)"""
        self.k_B = self.objSrouby.p_t/(2*Pi()) + self.objSrouby.mu_t * 0.9 * self.objSrouby.d_B0 / (2 * cos(self.objSrouby.alfa)) + \
            self.objSrouby.mu_n * 1.3 * self.objSrouby.d_B0

    def calcPhi_G(self):
        """(128)"""
        self.Phi_G = self.F_GI/(self.objTesneni.A_Gt *self.objTesneni.Q_smax)