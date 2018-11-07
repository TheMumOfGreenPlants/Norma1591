from math import pi, acos
from Priruba import *
from Sroub import *
from Tesneni import *
from Soucast import *
from Matice import *
import numpy

class Zatizeni(object):
    """description of class"""
    P_I = numpy.asarray([0,1])
    F_XI = numpy.asarray([0,0])
    F_YI = numpy.asarray([0,0])
    F_ZI = numpy.asarray([0,0])
    M_XI = numpy.asarray([0,0])
    M_YI = numpy.asarray([0,0])
    M_ZI = numpy.asarray([0,0])
    N_R = 1
        
    def setall(self, priruba1, priruba2, srouby, tesneni, matice ,podlozka1, podlozka2):
        self.objPriruba1 = priruba1
        self.objPriruba2 = priruba2
        self.objSrouby = srouby
        self.objTesneni = tesneni
        self.objMatice = matice
        self.objPodlozka1 = podlozka1
        self.objPodlozka2 = podlozka2


    def calcA_Q(self):
        """(90)"""
        self.objTesneni.calcb_Ge(self.F_G0)
        self.A_Q = (pi * self.objTesneni.d_Ge**2) /4     #d_Ge==d_Gi????

    def calcF_QI(self):
        """(91)"""
        self.calcA_Q()
        self.F_QI = self.A_Q*self.P_I

    def calcFM(self):
        """(92)(93)(94)(95)(96)"""
        self.F_AI = self.F_ZI
        self.F_LI = (self.F_XI**2 + self.F_YI**2)**(1/2)
        self.M_AI = (self.M_XI**2 + self.M_YI**2)**(1/2)
        self.M_TGI = self.M_ZI

        self.F_RI = [self.F_AI + (4 / self.objPriruba1.d_3e) * self.M_AI  ,\
                     self.F_AI - (4 / self.objPriruba1.d_3e) * self.M_AI]

    def calcdeltaU_TI(self):
        """(97)"""
        def part(soucast):
            return soucast.e * soucast.alfa * (soucast.T - soucast.T[0])
        self.deltaU_TI = self.objSrouby.l_B * self.objSrouby.alfa * (self.objSrouby.T - self.objSrouby.T[0]) \
            - part(self.objPriruba1 ) - part(self.objPriruba2) - part(self.objTesneni) - part(self.objPodlozka1) \
            - part(self.objPodlozka2)

    def conditionl_B(self):         
        """(98)"""
        sum = self.objPriruba1.e + self.objPriruba2.e + self.objPodlozka1.e + self.objPodlozka2.e + self.objTesneni.e                         
        return sum != self.objSrouby.l_B

    def calcY(self):
        """(99)(100)(101)(102)"""
        self.objTesneni.calcX_G()
        self.objPriruba1.calch_P(self.objTesneni)
        self.objPriruba2.calch_P(self.objTesneni)
        self.objPriruba1.calch_R()
        self.objPriruba2.calch_R()
        self.objPriruba1.calch_QGHL(self.objTesneni.d_Ge)
        self.objPriruba2.calch_QGHL(self.objTesneni.d_Ge)
        self.objSrouby.calcX_B()
        self.objPodlozka1.calcX_W(self.objPriruba1.d_5,self.objSrouby.d_B4,self.objSrouby.n_B)
        self.objPodlozka2.calcX_W(self.objPriruba2.d_5,self.objSrouby.d_B4,self.objSrouby.n_B)
        self.Y_BI = self.objPriruba1.Z_L * self.objPriruba1.h_L**2 / self.objPriruba1.E \
            + self.objPriruba2.Z_L * self.objPriruba2.h_L**2 / self.objPriruba2.E \
            + self.objSrouby.X_B / self.objSrouby.E + self.objPodlozka1.X_W / self.objPodlozka1.E \
            + self.objPodlozka2.X_W / self.objPodlozka2.E
        self.Y_GI = self.objPriruba1.Z_F * self.objPriruba1.h_G**2 / self.objPriruba1.E \
            + self.objPriruba2.Z_F * self.objPriruba2.h_G**2 / self.objPriruba2.E \
            + self.Y_BI + self.objTesneni.X_G/ self.objTesneni.E
        self.Y_QI = self.objPriruba1.Z_F * self.objPriruba1.h_G * (self.objPriruba1.h_H - self.objPriruba1.h_P + self.objPriruba1.h_Q)/self.objPriruba1.E \
            + self.objPriruba2.Z_F * self.objPriruba2.h_G * (self.objPriruba2.h_H - self.objPriruba2.h_P + self.objPriruba2.h_Q)/self.objPriruba2.E \
            + self.Y_BI
        self.Y_RI = self.objPriruba1.Z_F * self.objPriruba1.h_G * (self.objPriruba1.h_H + self.objPriruba1.h_R)/self.objPriruba1.E \
            + self.objPriruba2.Z_F * self.objPriruba2.h_G * (self.objPriruba2.h_H + self.objPriruba2.h_R)/self.objPriruba2.E \
            + self.Y_BI

    def calcF_G0min(self):
        """(103)"""
        self.F_G0min = self.objTesneni.A_Ge * self.objTesneni.Q_A

    def calcF_GImin(self):
        """(104)"""
        self.calcFM()
        self.objTesneni.calcd_Gt()
        self.objTesneni.calcA_Ge()
        self.calcF_G0min()
        F_GIminA = numpy.asarray([self.F_G0min])
        F_GIminB = numpy.asarray([self.F_G0min])
        for i in range(1,len(self.P_I)):
            F_GIminlocA = max(self.objTesneni.A_Ge * (self.objTesneni.Q_sminLI[i]),-((self.F_QI[i]) + (self.F_RI[0][i]) ),\
                        ((self.F_LI[i])/self.objTesneni.mu_G) + (2 * self.M_TGI[i]) / (self.objTesneni.mu_G * self.objTesneni.d_Gt) - (2 * self.M_AI[i]) / self.objTesneni.d_Gt)
            F_GIminlocB = max(self.objTesneni.A_Ge * self.objTesneni.Q_sminLI[i],-(self.F_QI[i] + self.F_RI[1][i]),\
                        self.F_LI[i]/self.objTesneni.mu_G + (2 * self.M_TGI[i]) / (self.objTesneni.mu_G * self.objTesneni.d_Gt) - (2 * self.M_AI[i]) / self.objTesneni.d_Gt)
            F_GIminA = numpy.append(F_GIminA, F_GIminlocA)
            F_GIminB = numpy.append(F_GIminB,F_GIminlocB)
        self.F_GImin = numpy.asarray([F_GIminA,F_GIminB])

    def calcdeltae_Gc(self):
        """(F.3)"""
        self.objTesneni.calcP_QR()
        self.deltae_Gc = self.objTesneni.K * self.Y_GI * self.objTesneni.deltae_Gc_test

    def calcF_GdeltaI(self):
        """(106) 1/2"""
        self.calcF_QI()
        self.calcF_GImin()
        self.calcY()
        self.calcdeltaU_TI()
        self.calcdeltae_Gc()
        arg = (self.F_GImin * self.Y_GI + self.F_QI * self.Y_QI + (self.F_RI * self.Y_RI - self.F_RI[0] * self.Y_RI[0]) \
                + self.deltaU_TI + self.deltae_Gc + (self.objTesneni.e_G - self.objTesneni.e_GA))
        arg = numpy.delete(arg,0,1)
        self.F_GdeltaI = arg/ self.Y_GI[0]

            ####tadyy!!!!!!!!!!!!


    def calcF_Gdelta(self):
        """(106) 2/2"""
        self.calcF_GdeltaI()
        self.F_Gdelta = max(self.F_GdeltaI)

    def calcF_G0req(self):
        """(107)"""
        self.F_G0 = 282018.6  # F_G0pocatecni - vlastni volba 
        self.F_G0 = self.objSrouby.A_B * self.objSrouby.f_B0 / 3 - self.objSrouby.F_R0
        self.F_G0req = 0
        while abs(self.F_G0req - self.F_G0) >= (self.F_G0req * 0.001):
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
        self.calcF_Gdelta()
        self.F_G0d = max(self.F_Gdelta , (2/3) * (1 - 10 / self.N_R) * self.F_B0max - self.F_R0)

    def calcF_GI(self):
        """(121)"""
        self.calcF_G0d()
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
        self.calcF_GI()
        self.Phi_G = self.F_GI/(self.objTesneni.A_Gt *self.objTesneni.Q_smax)