from math import *
from Soucast import *

class Tesneni(Soucast):
    """description of class"""
    Q_A = 100       # priloha G - neni pozadovana mira netesnosti   [MPa]
    d_G1 = 67       # teoreticky vnitrni prumer tesnici plochy      [mm]
    d_G2 = 120      # teoreticky vnejsi prumer tesnici plochy       [mm]
    e_G = 2         # tloustka tesneni v nezatizenem stavu          [mm]
    Q_smax = 30    # maximalni dovoleny tlak na tesneni            [MPa]
    Q_sminLI = numpy.asarray([8,8])    # minimalni povrchovy (utahovaci) tlak          [MPa]
                    # pusobici na tesneni , pozadovany pro tridu tesnosti L v podminkach provozu
    mu_G = 0.1
    # Zkouska P_QR
    T_PQR = numpy.asarray([20,160])
    Q_I = numpy.asarray([15,15])       # pocatecni napeti v tesneni                    [MPa]
    Q_R = numpy.asarray([13.7,5.3])       # zbytkove naapeti v tesneni                    [MPa]
    d_Gext = numpy.asarray([162.5,162.5])   # vnejsi prumer tesneni pouziteho pri zkousce   [mm]
    d_Gint = numpy.asarray([114.98,114.98])   # vnitrni prumer tesneni pouziteho pri zkousce  [mm]
    K = numpy.asarray([1500000,1500000])    # tuhost zk. zarizeni                           [N/mm]

    # Zkouska tloustky tesneni
    T_Gzk = numpy.asarray([20,160])
    Q_Gzk = numpy.asarray([[5,8,10,13,15,18,20,30,40,50,60,80],[5,2,3,4,4,5,6,7,10,13,17,20]])
    e_Gzk = numpy.asarray([[0.0119,0.0204,0.0296,0.0402,0.0519,0.0648,0.0789,0.1735,0.3031,0.4111,0.4958,0.6145],[0.0109,0.0815,0.1781,0.3119,0.4428,0.5525,0.6388,0.7117,0.9539,1.0841,1.1692,1.2312]])

    # Zkouska modulu pruznosti
    T_Ezk = numpy.asarray([20,160])
    Q_Ezk = numpy.asarray([[5,8,10,13,15,18,20,30,40,50,60,80],[5,2,3,4,4,5,6,7,10,13,17,20]])
    E_Ezk = numpy.asarray([[1120,1210,1194,1231,1284,1330,1396,1733,2032,2328,2664,3305],[1,94,90,104,132,176,220,118,213,315,413,390]])

    ##vypoctove parametry - uzivatel nemeni

    def tesneni_interp(self,T,Q,x,y,t):
        lin_k = numpy.zeros([len(t),2])
        for i in range(len(x)):
            lin_k[i] = (numpy.polyfit(x[i],y[i],1))
        X = (Q - lin_k[:,1]) / lin_k[:,0]

        X_fin = numpy.asarray([[Soucast.lin_interpolace(T,t,X[0])],[Soucast.lin_interpolace(T,t,X[1])]])
        return X_fin
     
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
        while numpy.all(numpy.absolute( b_Geold - self.b_Ge ) >= b_Geold * 0.001 ): 
            b_Geold = self.b_Ge
            self.calcb_Gi(obj1,obj2,F_G0)
            self.calcb_Ge(obj1,obj2,F_G0)
        self.calcX_G()
        obj1.calc645(self.d_Ge)
        obj2.calc645(self.d_Ge)

    def calce_G(self):
        """Interpolace tloustky tesneni dle krivky"""
        self.e = self.e_G - self.tesneni_interp(self.T[0],self.Q_G0,self.e_Gzk,self.Q_Gzk,self.T_Gzk)
                       
    def calcA_Ge(self):
        """(56)"""
        self.A_Ge = pi * self.d_Ge * self.b_Ge

    def calcQ_G0(self,F_G0):
        """(57)"""
        self.Q_G0 = F_G0 / self.A_Ge

    def calcE_G0(self,F_G0):
        """(58)"""
        self.calcQ_G0(F_G0)
        self.E_G0 = self.tesneni_interp(self.T[0],self.Q_G0,self.E_Ezk,self.Q_Ezk,self.T_Ezk)

    def getb_Gi(self):
        self.calcb_Gi()
        return self.b_Gi

    def calcX_G(self):
        """(63)"""
        self.X_G = (self.e_G/self.A_Gt)*(self.b_Gt + self.e_G/2) / (self.b_Ge + self.e_G/2)

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
        
