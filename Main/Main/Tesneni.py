from math import *
from PrirubaSKuzelovymKrkem import *
from Soucast import *

class Tesneni(Soucast):
    """description of class"""
    Q_A = 100         # priloha G - neni pozadovana mira netesnosti   [MPa]
    d_G1 = 67       # teoreticky vnitrni prumer tesnici plochy      [mm]
    d_G2 = 120      # teoreticky vnejsi prumer tesnici plochy       [mm]
    e_G = 2         # tloustka tesneni v nezatizenem stavu          [mm]
    #e_GA = e_G      # zjednoduseni!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Q_smax = 480    # maximalni dovoleny tlak na tesneni            [MPa]
    Q_sminLI = numpy.asarray([0,8])    # minimalni povrchovy (utahovaci) tlak          [MPa]
                    # pusobici na tesneni , pozadovany pro tridu tesnosti L v podminkach provozu
    mu_G = 0.1      #
    Q_I = 100       # pocatecni napeti v tesneni                    [MPa]
    Q_R = 100       # zbytkove naapeti v tesneni                    [MPa]
    d_Gext = 73.5   # vnejsi prumer tesneni pouziteho pri zkousce   [mm]
    d_Gint = 37.5   # vnitrni prumer tesneni pouziteho pri zkousce  [mm]
    K = 1500000     # tuhost zk. zarizeni                           [N/mm]
    typ = 1         # 1 - ploche kovove kruhove tesneni s pravouhlym prurezem; 2 - pro nekovova plocha tesneni
    ##vypoctove parametry - uzivatel nemeni

    def sete(self):
        self.e = self.e_G
        self.calce_GA()

    def calce_GA(self):
        self.e_GA = self.e

    def calcb_Gifirst(self): 
        """(64)"""
        self.b_Gi = self.b_Gt
        

    def calcb_Gt(self):
        """(51)"""
        self.b_Gt = (self.d_G2 - self.d_G1)/2

    def calcd_Gt(self):
        """(52)"""
        # teoreticky prumer tesneni
        self.d_Gt = (self.d_G1 + self.d_G2) / 2            

    def calcA_Gt(self):
        """(53)"""
        self.calcd_Gt()
        self.A_Gt = pi * self.d_Gt * self.b_Gt

    def calcb_Ge(self, F_G0):
        """(55)"""
        self.calcb_Gt()
        self.calcb_Gifirst()                
        self.b_Ge = min( self.b_Gi, self.b_Gt)
        
        self.calcb_Gi(F_G0)
        while numpy.all(numpy.absolute( self.b_Ge - self.b_Gi ) >= self.b_Ge * 0.001 ): 
            self.b_Ge = numpy.minimum( self.b_Gi, self.b_Gt)
            self.calcb_Gi(F_G0)
            self.b_Gi = numpy.minimum( self.b_Gi, self.b_Gt)
                       
    def calcA_Ge(self):
        """(56)"""
        self.A_Ge = pi * self.d_Ge * self.b_Ge

    def calcQ_G0(self):
        """(57)"""
        self.calcA_Ge()
        self.calcF_G0()
        self.Q_G0 = self.F_G0 / self.A_Ge

    def calcE_G0(self):
        """(58)"""
        self.E_G0 = self.E[0]

    def geth_G0(self, objPrirubaX):
        """(62)(61)(60)(59)"""
        if objPrirubaX == TocivaPrirubaSObrubou_Lemem: 
            h_G0 = (objPrirubaX.d_70 - self.d_Ge) / 2
        else:
            h_G0 = (objPrirubaX.d_3e - self.d_Ge) / 2
        return h_G0

    def getb_Gi(self):
        self.calcb_Gi()
        return self.b_Gi

    def calcX_G(self):
        """(63)"""
        self.calcA_Gt()
        self.X_G = (self.e_G/self.A_Gt)*(self.b_Gt + self.e_G/2) / (self.b_Ge + self.e_G/2)

    def calcb_Gi(self, F_G0):
        """(65)"""
        self.calcE_Gm()
        self.calcd_Ge()
        #self.b_Gi = self.b_Gt
        self.b_Gi = ( self.e / ( pi * self.d_Ge * self.E_Gm) / ( self.geth_G0(self.objPrvniPriruba) * self.objPrvniPriruba.Z_F / self.objPrvniPriruba.E[0]
                    + self.geth_G0(self.objDruhaPriruba) * self.objDruhaPriruba.Z_F / self.objDruhaPriruba.E[0] ) + 
                    ( F_G0 / (pi * self.d_Ge * self.Q_smax))**2)**(1/2)

    def calcE_Gm(self):
        """(66)(67)"""
        self.calcE_G0()
        if self.typ == 1:
            self.E_Gm = self.E_G0
        elif self.typ == 2:
            self.E_Gm = 0.5 * self.E_G0
        else:
            print('Zkontrolujte zvoleny typ tesneni!')
            sys.exit(int(0))

    def calcd_Ge(self):
        """tabulka 1 (68)"""
        self.d_Ge = self.d_G2 - self.b_Ge

    def calcF_G0min(self):
        """(103)"""
        self.calcA_Ge()
        self.F_G0min = self.A_Ge * self.Q_A

    def calcF_Gdelta(self):
        """(105) (106)"""
        self.F_Gdelta = 0 ## neuvazujeme jiny nez I=0 zatezny stav

    def calcF_G0req(self, F_G0):
        """(107)"""
        self.F_G0 = F_G0
        self.calcF_G0min()
        self.calcF_Gdelta()
        self.F_G0req = max(self.F_G0min,self.F_Gdelta)
        return(self.F_G0req)

    def setPriruby(self, objPriruba1, objPriruba2):
        self.objPrvniPriruba = objPriruba1
        self.objDruhaPriruba = objPriruba2

    def calcP_QR(self):
        """(F.1)(F.2)"""
        self.P_QR = self.Q_R / self.Q_I
        self.A_Gt_test = (pi / 4) * (self.d_Gext**2 - self.d_Gint**2)
        self.deltae_Gc_test = (self.A_Gt_test * self.Q_I * (1 - self.P_QR)) / self.K
        
    def calcPhi_G(self):
        """(128)"""
        self.calcPhi_B()
        self.calcF_G0max()
        self.F_G = numpy.insert(self.F_GI,[0],self.F_G0max,1) 
        self.F_G1 = numpy.insert(self.F_GI,[0],self.F_G0max1,1)
        self.F_GEXCEL = numpy.insert(self.F_GIEXCEL,[0],self.F_G0max1,1)
        self.Phi_G = self.F_G/(self.objTesneni.A_Gt *self.objTesneni.Q_smax)
        self.Phi_G1 = self.F_G1/(self.objTesneni.A_Gt *self.objTesneni.Q_smax)
        self.Phi_GEXCEL = self.F_GEXCEL/(self.objTesneni.A_Gt *self.objTesneni.Q_smax)



