from math import pi, acos
from Priruba import *
from Sroub import *
from Tesneni import *
from Soucast import *

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
    F_R0 = None
    Y_R0 = None
        
    def setall(self, priruba1, priruba2, srouby, tesneni):
        self.objPriruba1 = priruba1
        self.objPriruba2 = priruba2
        self.objSrouby = srouby
        self.objTesneni = tesneni

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

    def conditionl_B(self, e_Ft1, e_Ft2, e_L1, e_L2, e_Gt, e_W1, e_W2, l_B):
        """(98)"""
        sum = e_Ft1 + e_Ft2 + e_L1 + e_L2 + e_Gt + e_W1 + e_W2
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
        if Y_R0 == None:
            self.Y_R0 = self.Y_RI
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
