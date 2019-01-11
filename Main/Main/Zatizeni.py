from math import pi, acos
from Priruba import *
from Sroub import *
from Tesneni import *
from Soucast import *
from Matice import *
import numpy
import sys

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
        self.A_Q = (pi * self.objTesneni.d_G1**2) /4     # v norme je d_Ge, ale pry se nahrazuje !!!d_G1????!!!!

    def calcF_QI(self):
        """(91)"""
        self.calcA_Q()
        #F_QIA = self.A_Q[0] * self.P_I
        #F_QIB = self.A_Q[1] * self.P_I
        #self.F_QI = numpy.asarray([F_QIA,F_QIB])
        self.F_QI = self.A_Q * self.P_I

    def calcFM(self):
        """(92)(93)(94)(95)(96)"""
        self.F_AI = self.F_ZI
        self.F_LI = (self.F_XI**2 + self.F_YI**2)**(1/2)
        self.M_AI = (self.M_XI**2 + self.M_YI**2)**(1/2)
        self.M_TGI = self.M_ZI

        self.F_RI = numpy.asarray([self.F_AI + (4 / self.objPriruba1.d_3e) * self.M_AI  ,\
                     self.F_AI - (4 / self.objPriruba1.d_3e) * self.M_AI])

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
        self.objPriruba1.calch_P(self.objTesneni)
        self.objPriruba2.calch_P(self.objTesneni)
        self.objPriruba1.calch_R()
        self.objPriruba2.calch_R()
        self.objPriruba1.calch_QGHL(self.objTesneni.d_Ge)
        self.objPriruba2.calch_QGHL(self.objTesneni.d_Ge)
        self.objSrouby.calcX_B()                                #ok
        self.objTesneni.calcX_G()
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
        self.objTesneni.calcd_Gt()
        self.objTesneni.calcA_Ge()
        self.calcF_G0min()
        self.F_GImin = numpy.maximum(numpy.asarray([self.objTesneni.A_Ge[0] * self.objTesneni.Q_sminLI,self.objTesneni.A_Ge[1] * self.objTesneni.Q_sminLI]),\
              - (self.F_QI + self.F_RI),numpy.asarray([self.F_LI/self.objTesneni.mu_G + (2 * self.M_TGI) / (self.objTesneni.mu_G * self.objTesneni.d_Gt) \
              - (2 * self.M_AI) / self.objTesneni.d_Gt,self.F_LI/self.objTesneni.mu_G + (2 * self.M_TGI) / (self.objTesneni.mu_G * self.objTesneni.d_Gt) \
              - (2 * self.M_AI) / self.objTesneni.d_Gt]))


    def calcdeltae_Gc(self):
        """(F.3)"""
        self.objTesneni.calcP_QR()
        self.deltae_Gc = self.objTesneni.K * self.Y_GI * self.objTesneni.deltae_Gc_test

    def calcF_GdeltaI(self):
        """(106) 1/2"""
        self.calcF_QI()         #ok
        self.calcF_GImin()      #ok
        self.calcY()            #ok
        self.calcdeltaU_TI()
        self.calcdeltae_Gc()

        arg = self.F_GImin[:,1:] * self.Y_GI[1:] + self.F_QI[1:] * self.Y_QI[1:] + self.deltaU_TI[1:] + + self.deltae_Gc[1:] + (self.objTesneni.e_G - self.objTesneni.e_GA)
        self.F_GdeltaI = arg/ self.Y_GI[0]

    def calcF_Gdelta(self):
        """(106) 2/2"""
        self.calcF_GdeltaI()
        self.F_Gdelta = numpy.vstack((numpy.max(self.F_GdeltaI,1)))

    def setF_G0(self):
        """(1)(54)"""
        self.objSrouby.calcEps()
        if self.objSrouby.zatizeni_sroubu:
            self.calcFM()
            self.F_G0 = numpy.vstack((self.objSrouby.F_B0spec * (1-self.objSrouby.Eps_minus) - self.F_RI[:,0]))     #(1)
        else:
            self.F_G0 = numpy.vstack((self.objSrouby.A_B * self.objSrouby.f_B0 / 3 - self.F_RI[:,0]))               #(54


    def calcF_G0req(self):
        """(107)"""
        self.calcFM()
        self.setF_G0()
        self.F_G0req = numpy.asarray([numpy.asarray([0]),numpy.asarray([0])])
        self.calcF_Gdelta()
        self.F_G0req = numpy.vstack((numpy.maximum(numpy.concatenate(self.F_G0min),numpy.concatenate(self.F_Gdelta))))
        A=(self.F_G0req > self.F_G0)
        B=self.objSrouby.zatizeni_sroubu
        C=numpy.all(A.all and B)
        if (self.objSrouby.zatizeni_sroubu and (self.F_G0req > self.F_G0).all):
            print('Pro splneni kriterii tesnosti se musi zvysit hodnota F_B0spec a vypocet spustit znovu!')
            sys.exit(int(0))
        else:
            while numpy.all(numpy.absolute(self.F_G0req - self.F_G0) >= (self.F_G0req * 0.001)):
                self.calcF_Gdelta()
                self.F_G0req = numpy.vstack((numpy.maximum(numpy.concatenate(self.F_G0min),numpy.concatenate(self.F_Gdelta))))
                self.F_G0 = self.F_G0req
                self.calcF_Gdelta()
                self.F_G0req = numpy.vstack((numpy.maximum(numpy.concatenate(self.F_G0min),numpy.concatenate(self.F_Gdelta))))


    def calcF_B0req(self):
        """(108)"""
        self.calcF_G0req()
        self.F_B0req = self.F_G0req + numpy.vstack((self.F_RI[:,0]))

    def calcF_B0minmax(self):
        """(112)(113)"""
        self.calcF_B0av()
        self.F_B0min = self.F_B0av * (1 - self.objSrouby.Eps_minus)
        self.F_B0max = self.F_B0av * (1 - self.objSrouby.Eps_plus)

    def calcF_B0nom(self):
        """(114)"""
        self.calcF_B0req()  
        self.calcF_B0minmax()
        if self.objSrouby.kontrola:
            self.calcF_B0nom1()
        else:
            self.calcF_B0nom2()
        self.P = self.F_B0min >= self.F_B0req

    def calcF_B0nom1(self):
        """(115)"""
        #metoda utahovani sroubu S kontrolou zatizeni sroubu
        self.F_B0nom = self.F_B0req/(1-self.objSrouby.Eps_minus)
        
    def calcF_B0av(self):
        """(B.3)"""
        self.F_B0av = min(self.objSrouby.A_B * self.objSrouby.f_B0, self.objSrouby.n_B *200000 )
        
    def calcF_B0nom2(self):
        """(116)"""
        #metoda utahovani sroubu BEZ kontroly zatizeni sroubu
        self.objSrouby.Eps1_minus = 0.5
        self.objSrouby.calcEps()
        self.P_A2 = self.F_B0av >= self.F_B0req/(1-self.objSrouby.Eps_minus)
        self.F_B0nom = self.F_B0av

    def calcF_B0max(self):
        """(117)"""
        self.calcF_B0nom()
        self.F_B0max = self.objSrouby.F_B0nom * ( 1 + self.objSrouby.Eps_plus )

    def calcF_G0max(self):
        """(118)"""
        self.F_G0max = self.F_B0max - self.F_R0

    def calcF_G0d(self):
        """(2)(119)"""
        if self.objSrouby.zatizeni_sroubu:
            self.F_G0d = numpy.maximum(self.F_B0min - self.F_RI[:,0] , (2/3) * (1 - 10 / self.N_R) * self.F_B0max - self.F_RI[:,0])         #(2) pri vypoctu s predem znamym zatizenim sroubu F_B0spec
        else:
            self.F_G0d = numpy.maximum(self.F_Gdelta , numpy.vstack((2/3) * (1 - 10 / self.N_R) * self.F_B0max - self.F_RI[:,0]))                         #(119) v ostatnich pripadech
       
    def calcF_GI(self):
        """(121)"""
        self.calcF_G0d()
        A=self.F_G0d * numpy.vstack((self.Y_GI[:,0]))
        B=self.F_QI[1:] * self.Y_QI[:,1:]
        C=self.F_RI[:,1:] * self.Y_RI[:,1:]
        D=numpy.vstack((self.F_RI[:,0] * self.Y_RI[:,0]))
        E=self.deltaU_TI[1:]
        F=self.deltae_Gc[:,1:]
        G=self.objTesneni.e_G - self.objTesneni.e_GA
        H=numpy.vstack((self.Y_GI[:,0]))
        self.F_GI = (self.F_G0d * numpy.vstack((self.Y_GI[:,0])) - (self.F_QI[1:] * self.Y_QI[:,1:] + (self.F_RI[:,1:] * self.Y_RI[:,1:] - numpy.vstack((self.F_RI[:,0] * self.Y_RI[:,0]))) \
                + self.deltaU_TI[1:] ) - self.deltae_Gc[:,1:] - (self.objTesneni.e_G - self.objTesneni.e_GA)) / numpy.vstack((self.Y_GI[:,1:]))
        F_G0nom = (self.F_B0nom - numpy.vstack((self.F_RI[:,0])))
        self.F_GIEXCEL = (F_G0nom*numpy.vstack((self.Y_GI[:,0]))*self.objTesneni.P_QR-self.Y_RI[:,1:]*self.F_RI[:,1:]+numpy.vstack((self.F_RI[:,0]*self.Y_RI[:,0]))-self.F_QI[1:]*self.Y_QI[:,1:] \
            -self.deltaU_TI[1:])/numpy.vstack((self.Y_GI[:,0]))

    def calcF_BI(self):
        """(122)"""
        self.calcF_GI()
        self.F_BI=(self.F_GI + (self.F_QI[1:] + self.F_RI[:,1:]))
        self.F_BIEXCEL = (self.F_GIEXCEL + (self.F_QI[1:]  + self.F_RI[:,1:] ))

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

    def calcPhi_B(self):
        """(123)"""
        self.calcM_tBnom()
        self.calcF_BI()
        self.calcc_A()
        self.calcc_B()
        self.M_tBnomEXCEL[0]=989.29
        self.M_tBnomEXCEL[1]=989.29
        self.F_BIEXCEL = numpy.insert(self.F_BIEXCEL,[0],self.F_B0nom,1)
        
        self.Phi_B = (1 / (self.objSrouby.f_B0 * self.c_B)) * ((self.F_BI/self.objSrouby.A_B)**2 + \
           3*(self.c_AI * self.M_tBnom*1000 / self.objSrouby.l_B)**2)**(1/2)
        self.Phi_BEXCEL = ((1 / (self.objSrouby.f_B0 * self.c_B))) * ((self.F_BIEXCEL/self.objSrouby.A_B)**2 + \
           3*(self.c_AI * self.M_tBnomEXCEL / (pi*self.objSrouby.d_Bs**3/16))**2)**(1/2)
        a=1

    def calcM_tBnom(self):
        """(B.9)"""
        self.objPodlozka1.calcX_W(self.objPriruba1.d_5,self.objSrouby.d_B4,self.objSrouby.n_B)
        self.objSrouby.calck_B(self.objPodlozka1.d_n)
        self.calcF_B0nom()
        self.M_tBnom = ((0.159*self.objSrouby.p_t+0.577*self.objSrouby.mu_t*self.objSrouby.d_B0*0.9) * self.F_B0nom / self.objSrouby.n_B)/1000  #/1000->Nm
        self.M_tBnomEXCEL = numpy.concatenate((0.159*self.objSrouby.p_t+0.519*self.objSrouby.mu_t*self.objSrouby.d_B0) * self.F_B0nom )/self.objTesneni.E
        self.M_tnom = self.objSrouby.k_B * self.F_B0nom / self.objSrouby.n_B
        self.M_tnomEXCEL = self.objSrouby.k_B * self.F_B0nom / 2335

    def calcPhi_G(self):
        """(128)"""
        self.calcF_GI()
        self.Phi_G = self.F_GI/(self.objTesneni.A_Gt *self.objTesneni.Q_smax)

    def calcPreload(self):                                                                                 # vypocet predpeti ve sroubu
        self.calcF_B0nom()                                                                                   # k vypoctu potrebujeme znat F_B0nom
        self.Preload = self.F_B0nom / self.objSrouby.A_B                                                                      # predpeti ve sroubu                                            [MPa]

