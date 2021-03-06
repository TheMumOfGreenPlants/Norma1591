from math import *
from Soucast import *

class Tesneni(Soucast):
    """description of class"""
    Q_A = 0       # priloha G - neni pozadovana mira netesnosti   [MPa]
    d_G1 = 0       # teoreticky vnitrni prumer tesnici plochy      [mm]
    d_G2 = 0      # teoreticky vnejsi prumer tesnici plochy       [mm]
    e_G = 0         # tloustka tesneni v nezatizenem stavu          [mm]
    Q_smax = 0    # maximalni dovoleny tlak na tesneni            [MPa]
    Q_sminLI = numpy.asarray([0,0])    # minimalni povrchovy (utahovaci) tlak          [MPa]
                    # pusobici na tesneni , pozadovany pro tridu tesnosti L v podminkach provozu
    mu_G = 0
    # Zkouska P_QR
    T_PQR = numpy.asarray([0,0])
    Q_I = numpy.asarray([0,0])       # pocatecni napeti v tesneni                    [MPa]
    Q_R = numpy.asarray([0,0])       # zbytkove naapeti v tesneni                    [MPa]
    d_Gext = numpy.asarray([0,0])   # vnejsi prumer tesneni pouziteho pri zkousce   [mm]
    d_Gint = numpy.asarray([0,0])   # vnitrni prumer tesneni pouziteho pri zkousce  [mm]
    K = numpy.asarray([0,0])    # tuhost zk. zarizeni                           [N/mm]

    # Zkouska tloustky tesneni
    T_Gzk = numpy.asarray([0,0])
    Q_Gzk = numpy.asarray([[0],[0]])
    e_Gzk = numpy.asarray([[0],[0]])

    # Zkouska modulu pruznosti
    T_Ezk = numpy.asarray([0,0])
    Q_Ezk = numpy.asarray([[0],[0]])
    E_Ezk = numpy.asarray([[0],[0]])

    ##vypoctove parametry - uzivatel nemeni

    def sete(self):
        self.e = self.e_G
        self.e_GA = self.e_G

    def calc642(self):
        """(51)(52)(53)"""
        self.b_Gt = (self.d_G2 - self.d_G1) / 2
        self.d_Gt = (self.d_G1 + self.d_G2) / 2         # teoreticky prumer tesneni
        self.A_Gt = pi * self.d_Gt * self.b_Gt
        self.calcP_QR()

    def calc643(self,obj1,obj2,F_G0):
        self.calcb_Gifirst()
        self.calcb_Ge(obj1,obj2,F_G0)

    def calc_mat_param(self):
        self.E = numpy.zeros(len(self.T))
        for t_i,t in enumerate(self.T):
            self.E[t_i] = Soucast.double_interp(t,self.Q_G0,self.E_Ezk,self.Q_Ezk,self.T_Ezk)
        self.alfa = Soucast.calc_mat_parameter(self.T,self.T_azk,self.alfa_zk)

    def calcb_Ge(self,obj1,obj2,F_G0):
        """(55)"""
        self.b_Ge = numpy.minimum( self.b_Gi, self.b_Gt)
        self.calcd_Ge()
        self.calcA_Ge()
        self.calcE_G0(F_G0)
        obj1.calch_G0(self.d_Ge)
        obj2.calch_G0(self.d_Ge)
        self.calce_G()
        self.calcE_Gm(F_G0)


    def iteraceb(self,obj1,obj2,F_G0):
        self.calcb_Gi(obj1,obj2,F_G0)
        b_Geold = self.b_Ge
        self.calcb_Ge(obj1,obj2,F_G0)
        while (numpy.absolute( b_Geold - self.b_Ge ) >= b_Geold * 0.001 ): 
            b_Geold = self.b_Ge
            self.calcb_Gi(obj1,obj2,F_G0)
            self.calcb_Ge(obj1,obj2,F_G0)
        self.calcX_G()
        obj1.calc645(self.d_Ge)
        obj2.calc645(self.d_Ge)
        self.calc_mat_param()

    def calce_G(self):
        """Interpolace tloustky tesneni dle krivky"""
        self.e = self.e_G - Soucast.double_interp(self.T[0],self.Q_G0,self.e_Gzk,self.Q_Gzk,self.T_Gzk)
                       
    def calcA_Ge(self):
        """(56)"""
        self.A_Ge = pi * self.d_Ge * self.b_Ge

    def calcQ_G0(self,F_G0):
        """(57)"""
        self.Q_G0 = F_G0 / self.A_Ge

    def calcE_G0(self,F_G0):
        """(58)"""
        self.calcQ_G0(F_G0)
        self.E_G0 = Soucast.double_interp(self.T[0],self.Q_G0,self.E_Ezk,self.Q_Ezk,self.T_Ezk)

    def getb_Gi(self):
        self.calcb_Gi()
        return self.b_Gi

    def calcX_G(self):
        """(63)"""
        self.X_G = (self.e/self.A_Gt)*(self.b_Gt + self.e/2) / (self.b_Ge + self.e/2)

    def calcF_G0min(self):
        """(103)"""
        self.calcA_Ge()
        self.F_G0min = self.A_Ge * self.Q_A

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
        
