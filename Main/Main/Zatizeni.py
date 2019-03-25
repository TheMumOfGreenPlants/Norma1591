from math import pi, acos
import numpy
import sys
from Soucast import *

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

    def setTypVypoctu(self, typ):
        self.typ = typ

    def setall(self, priruba1, priruba2, srouby, tesneni, matice ,podlozka1, podlozka2):
        self.objPriruba1 = priruba1
        self.objPriruba2 = priruba2
        self.objSrouby = srouby
        self.objTesneni = tesneni
        self.objMatice = matice
        self.objPodlozka1 = podlozka1
        self.objPodlozka2 = podlozka2

    def calc722(self):
        self.calcA_Q()
        self.calcF_QI()
        self.calcdeltaU_TI()

    def calc74(self):
        self.calcdeltae_Gc()
        self.calcF_G0min()
        self.calcF_GImin()

    def calc75(self):
        self.calcF_Gdelta()

    def calc7(self):
        self.calc722()
        self.calc73()
        self.calc74()
        self.calc75()

    def calc752(self):
        self.calcF_B0minmax()
        self.calcF_B0nom()
        self.P1 = self.F_B0min >= self.F_B0req
        self.calcF_B0max()
        self.calcF_G0max()

    def calc76(self):
        self.calcF_G0d()
        self.calcF_GI()
        self.calcF_BI()

    def calc8(self):
        self.calcM_tBnom()
        self.F_B = numpy.insert(self.F_BI,[0],self.F_B0max,1)
        self.F_B1 = numpy.insert(self.F_BI,[0],self.F_B0max1,1)
        self.F_BEXCEL = numpy.insert(self.F_BIEXCEL,[0],self.F_B0nom,1)

    def iteraceF(self):
        self.calcF_G0req()
        if (self.typ == "KontrolaSroubu" and (self.F_G0req > self.F_G0).any):
            print('Pro splneni kriterii tesnosti se musi zvysit hodnota F_B0spec a vypocet spustit znovu!')
            sys.exit(int(0))
        while numpy.all(numpy.absolute(self.F_G0req - self.F_G0) >= (self.F_G0req * 0.001)):
            self.F_G0 = self.F_G0req
            self.objTesneni.iteraceb(self.objPriruba1,self.objPriruba2,self.F_G0)
            self.calcF_G0req()
        self.calcF_B0req()

    def calcA_Q(self):
        """(90)"""
        A_Qnorm = (pi * self.objTesneni.d_Ge**2) /4
        self.A_Q = (pi * self.objTesneni.d_G1**2) /4     # v norme je d_Ge, ale pry se nahrazuje !!!d_G1????!!!!
        
    def calcF_QI(self):
        """(91)"""
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

    def calc73(self):
        """(99)(100)(101)(102)"""
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
        self.F_GImin = numpy.maximum(numpy.maximum(self.objTesneni.A_Ge* self.objTesneni.Q_sminLI,- (self.F_QI + self.F_RI)),\
            self.F_LI/self.objTesneni.mu_G + (2 * self.M_TGI) / (self.objTesneni.mu_G * self.objTesneni.d_Gt) - (2 * self.M_AI) / self.objTesneni.d_Gt)[:,1:]

    def calcdeltae_Gc(self):
        """(F.3)"""
        deltae_Gc = self.objTesneni.K * self.Y_GI * self.objTesneni.deltae_Gc_test
        self.deltae_Gc = numpy.zeros([len(self.objTesneni.T),2])
        for i in range(len(self.objTesneni.T)):
            for j in range(2):
                self.deltae_Gc[j,i] = Soucast.lin_interpolace(self.objTesneni.T[i],self.objTesneni.T_PQR, deltae_Gc[j])

    def calcF_Gdelta(self):
        """(106)"""
        ### po doprogramovani e_GA odstranit:
        self.objTesneni.e_GA = self.objTesneni.e
        ###

        arg = self.F_GImin * self.Y_GI[:,1:] + self.F_QI[1:] * self.Y_QI[:,1:] +\
            self.deltaU_TI[1:] + self.deltae_Gc[:,1:] + (self.objTesneni.e - self.objTesneni.e_GA)
        self.F_Gdelta = numpy.asarray([[max(arg[0,:] / self.Y_GI[0,0])],\
                                      [max(arg[1,:] / self.Y_GI[1,0])]])

    def setF_G0(self):
        """(1)(54)"""
        self.objSrouby.calcEps()
        self.calcFM()
        if self.typ == "KontrolaSroubu":
            self.F_G0 = numpy.vstack((self.objSrouby.F_B0spec * (1-self.objSrouby.Eps_minus) - self.F_RI[:,0]))     #(1)
        elif self.typ == "MiraNetesnosti":
            self.F_G0 = numpy.vstack((self.objSrouby.A_B * self.objSrouby.f_B0 / 3 - self.F_RI[:,0]))               #(54)


    def calcF_G0req(self):
        """(107)(109)"""
        self.F_G0req = numpy.vstack(([[max(self.F_G0min[0],self.F_Gdelta[0])],\
                                    [max(self.F_G0min[1],self.F_Gdelta[1])]]))

    def calcF_B0req(self):
        """(108)"""
        self.F_B0req = self.F_G0req + numpy.vstack((self.F_RI[:,0]))

    def calcF_B0minmax(self):
        """(112)(113)"""
        self.calcF_B0av()
        self.F_B0min = self.F_B0av * (1 - self.objSrouby.Eps_minus)
        self.F_B0max = self.F_B0av * (1 - self.objSrouby.Eps_plus)


    def calcF_B0nom(self):
        """(114)"""
        if self.objSrouby.kontrola:
            self.calcF_B0nom1()
        else:
            self.calcF_B0nom2()


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
        self.F_B0max1 = self.F_B0nom * ( 1 + self.objSrouby.Eps_plus )

    def calcF_G0max(self):
        """(118)"""
        self.F_G0max1 = self.F_B0max1 - numpy.vstack((self.F_RI[:,0]))

    def calcF_G0d(self):
        """(2)(119)"""
        if self.typ == "KontrolaSroubu":
            self.F_G0d = numpy.maximum(self.F_B0min - self.F_RI[:,0] , (2/3) * (1 - 10 / self.N_R) * self.F_B0max - self.F_RI[:,0])         # (2) pri vypoctu s predem znamym zatizenim sroubu F_B0spec
        elif self.typ == "MiraNetesnosti":
            self.F_G0d = numpy.maximum(self.F_Gdelta , numpy.vstack((2/3) * (1 - 10 / self.N_R) * self.F_B0max - self.F_RI[:,0]))           # (119) v ostatnich pripadech
            self.F_G0dopt = numpy.maximum(self.F_B0min - self.F_RI[:,0] , (2/3) * (1 - 10 / self.N_R) * self.F_B0max - self.F_RI[:,0])
       
    def calcF_GI(self):
        """(121)"""
        self.F_GI = (self.F_G0d * numpy.vstack((self.Y_GI[:,0])) - (self.F_QI[1:] * self.Y_QI[:,1:] + (self.F_RI[:,1:] * self.Y_RI[:,1:] - numpy.vstack((self.F_RI[:,0] * self.Y_RI[:,0]))) \
                + self.deltaU_TI[1:] ) - self.deltae_Gc[:,1:] - (self.objTesneni.e - self.objTesneni.e_GA)) / numpy.vstack((self.Y_GI[:,1:]))
        self.F_G0nomEXCEL = (self.F_B0nom - numpy.vstack((self.F_RI[:,0])))
        self.F_GIEXCEL = (self.F_G0nomEXCEL*numpy.vstack((self.Y_GI[:,0]))*self.objTesneni.P_QR-self.Y_RI[:,1:]*self.F_RI[:,1:] + numpy.vstack((self.F_RI[:,0]*self.Y_RI[:,0]))-self.F_QI[1:]*self.Y_QI[:,1:] \
            -self.deltaU_TI[1:])/numpy.vstack((self.Y_GI[:,0]))

    def calcF_BI(self):
        """(122)"""
        self.F_BI=(self.F_GI + (self.F_QI[1:] + self.F_RI[:,1:]))
        self.F_BIEXCEL = (self.F_GIEXCEL + (self.F_QI[1:]  + self.F_RI[:,1:] ))



    def calcM_tBnom(self):
        """(B.9)"""
        self.M_tBnom = ((0.159*self.objSrouby.p_t+0.577*self.objSrouby.mu_t*self.objSrouby.d_B0*0.9) * self.F_B0nom / self.objSrouby.n_B)/1000  #/1000->Nm
        self.M_tBnomEXCEL = numpy.concatenate((0.159*self.objSrouby.p_t+0.519*self.objSrouby.mu_t*self.objSrouby.d_B0) * self.F_B0nom ) / 2335
        self.objSrouby.calck_B()
        self.M_tnom = self.objSrouby.k_B * self.F_B0nom / self.objSrouby.n_B / 1000
        self.M_tnomEXCEL = self.objSrouby.k_B * self.F_B0nom / 2335

    def calcTheta_F(self,p):
        """(C.1)(C.2)(C.3)(C.4)(C.5)(C.6)(C.7)(C.8)(C.9)(C.10)"""
        F_B0min = self.F_B0nom*(1-self.objSrouby.Eps_minus)
        F_B0max = self.F_B0nom*(1+self.objSrouby.Eps_plus)
        F_G0min = F_B0min - numpy.vstack((self.F_RI[:,0]))
        F_G0max = F_B0max - numpy.vstack((self.F_RI[:,0]))
        F_GImin = (F_G0min * numpy.vstack((self.Y_GI[:,0])) - (self.F_QI[1:] * self.Y_QI[:,1:] + (self.F_RI[:,1:] * self.Y_RI[:,1:] - numpy.vstack((self.F_RI[:,0] * self.Y_RI[:,0]))) \
                + self.deltaU_TI[1:] ) - self.deltae_Gc[:,1:]) / numpy.vstack((self.Y_GI[:,1:]))
        F_GImax = (F_G0max * numpy.vstack((self.Y_GI[:,0])) - (self.F_QI[1:] * self.Y_QI[:,1:] + (self.F_RI[:,1:] * self.Y_RI[:,1:] - numpy.vstack((self.F_RI[:,0] * self.Y_RI[:,0]))) \
                + self.deltaU_TI[1:] ) - self.deltae_Gc[:,1:]) / numpy.vstack((self.Y_GI[:,1:]))
        F_BImin = F_GImin + (self.F_QI[1:] +self.F_RI[1:])
        F_BImax = F_GImax + (self.F_QI[1:] +self.F_RI[1:])

        F_Gmin = numpy.insert(F_GImin,[0],self.F_G0min,1)
        F_Gmax = numpy.insert(F_GImax,[0],self.F_G0max,1)
        F_Bmin = numpy.insert(F_GImin,[0],self.F_B0min,1)
        F_Bmax = numpy.insert(F_GImax,[0],self.F_B0max,1)

        Theta_Fmin = (p.Z_F/p.E)*(F_Gmin * p.h_G + self.F_QI * (p.h_H - p.h_P + p.h_Q) + self.F_RI * (p.h_H - p.h_R))*(360/(2*pi))
        Theta_Fmax = (p.Z_F/p.E)*(F_Gmax * p.h_G + self.F_QI * (p.h_H - p.h_P + p.h_Q) + self.F_RI * (p.h_H - p.h_R))*(360/(2*pi))

        Theta_Lmin = (p.Z_L/p.E) * F_Bmin * p.h_L
        Theta_Lmax = (p.Z_L/p.E) * F_Bmax * p.h_L

        self.Theta_F = (p.Z_F/p.E)*(self.F_G * p.h_G + self.F_QI * (p.h_H - p.h_P + p.h_Q) + self.F_RI * (p.h_H - p.h_R))*(360/(2*pi))
        self.Theta_F1 = (p.Z_F/p.E)*( (self.F_G1 * p.h_G) + self.F_QI * (p.h_H - p.h_P + p.h_Q) + self.F_RI * (p.h_H - p.h_R))*(360/(2*pi))

    def calcPhi_B(self):
        """(123)"""
        self.calcc_A()
        self.calcc_B()
        
        self.Phi_B = (1 / (self.objSrouby.f_B0 * self.c_B)) * ((self.F_B/self.objSrouby.A_B)**2 + \
           3*(self.c_AI * self.M_tBnom*1000 / self.objSrouby.l_B)**2)**(1/2)
        self.Phi_B1 = (1 / (self.objSrouby.f_B0 * self.c_B)) * ((self.F_B1/self.objSrouby.A_B)**2 + \
           3*(self.c_AI * self.M_tBnom*1000 / (pi*self.objSrouby.d_Bs**3/16))**2)**(1/2)
        self.Phi_B2 = (1 / (self.objSrouby.f_B0 * self.c_B)) * ((self.F_B/self.objSrouby.A_B)**2 + \
           3*(self.c_AI * self.M_tBnom*1000 / (pi*self.objSrouby.d_Bs**3/16))**2)**(1/2)
        self.Phi_BEXCEL = ((1 / (self.objSrouby.f_B0 * self.c_B))) * ((self.F_BEXCEL/self.objSrouby.A_B)**2 + \
           3*(self.c_AI * self.M_tnomEXCEL / (pi*self.objSrouby.d_Bs**3/16))**2)**(1/2)
        x=1

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

    def calcPhi_G(self):
        """(128)"""
        self.F_GEXCEL = numpy.insert(self.F_GIEXCEL,[0],self.F_G0max1,1)
        self.Phi_G = self.F_G/(self.objTesneni.A_Gt *self.objTesneni.Q_smax)
        self.Phi_G1 = self.F_G1/(self.objTesneni.A_Gt *self.objTesneni.Q_smax)
        self.Phi_GEXCEL = self.F_GEXCEL/(self.objTesneni.A_Gt *self.objTesneni.Q_smax)