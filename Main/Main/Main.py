from Priruba import *
from ObecnaPriruba import *
from ZaslepovaciPriruba import *
from TocivaPrirubaSObrubou_Lemem import *



from Podlozka import *
from Sroub import *
from Matice import *
from Tesneni import *
from TesneniTyp1 import *
from TesneniTyp2 import *
from TesneniTyp3 import *
from TesneniTyp4 import *
from Zatizeni import *
from math import pi, acos
from sys import stdin
import sys
import copy

def main():
    # Volby uzivatele:
    def VolbaPriruby(typ):
        return {
            1 : ObecnaPriruba(1,1,1),
            2 : ObecnaPriruba(1,1,2),
            3 : ObecnaPriruba(0,2,1),
            4 : ObecnaPriruba(0,2,2),
            5 : ObecnaPriruba(0,3,1),
            6 : ObecnaPriruba(0,3,2),
            7 : ObecnaPriruba(0,1,1),
            8 : ZaslepovaciPriruba(),
            9 : TocivaPrirubaSObrubou_Lemem(0,1),
            10: TocivaPrirubaSObrubou_Lemem(0,0),
            11: ObecnaPriruba(1,1,1),
            12: TocivaPrirubaSObrubou_Lemem(1,0),
            }[typ]
    objSrouby = Sroub()
    def VolbaTesneni(typ):
        A = {
            1 : TesneniTyp1,
            2 : TesneniTyp2,
            3 : TesneniTyp3,
            4 : TesneniTyp4,
            }[typ]
        return A()
    objMatice = Matice()
    objZatizeni = Zatizeni()
    # Volba typu vypoctu
    def ZvolTypVypoctu(typ):
        return {
            1 : "MiraNetesnosti",
            2 : "KontrolaSroubu",
            }[typ]


    #################################################################
    ###################    NASTAVENI PARAMETRU   ####################
    #################################################################

    TypVypoctu = ZvolTypVypoctu(1)

    objPrvniPodlozka = Podlozka(0)
    objDruhaPodlozka = Podlozka(0)

    ##### objPrvniPriruba:  ###### 
    objPrvniPriruba = VolbaPriruby(7)
    # Pro vsechny typy prirub:
    objPrvniPriruba.d_0 = 110
    objPrvniPriruba.d_3 = 170
    objPrvniPriruba.d_4 = 200
    objPrvniPriruba.d_5 = 19
    objPrvniPriruba.e_Fb = 16.5
    objPrvniPriruba.e_Ft = objPrvniPriruba.e_Fb 
    objPrvniPriruba.e_F = objPrvniPriruba.e_Fb       # vypocet dle 2 * A_F /(d_4 - d_0)
    objPrvniPriruba.Fi_S = 0        # 0 - pro valec, natoceni pripojne skorepiny                    [rad]
    objPrvniPriruba.d_S = 112.1        # #POUZE priruba bez krku; stredni prumer skorepiny (prumer v miste spoje s prirubou)        [mm]
    objPrvniPriruba.f_F = 205         # dovolene namahani priruby     [MPa]
    objPrvniPriruba.f_S = 205         # dovolene namahani skorepiny   [MPa]

    # Pro Obecnou Prirubu:
    objPrvniPriruba.e_1 = 2.1     # nejmensi tloustka steny na tenkem konci krku   [mm]
    

    objPrvniPriruba.T_Ezk = numpy.asarray([20,100])
    objPrvniPriruba.E_zk = numpy.asarray([200000,200000])
    objPrvniPriruba.T_azk = numpy.asarray([20,100])
    objPrvniPriruba.alfa_zk = numpy.asarray([12.6e-6,12.6e-6])

    
    ##### objDruhaPriruba:  ###### 
    objDruhaPriruba = copy.deepcopy(objPrvniPriruba)

    ##### objSrouby:  ###### 
    objSrouby.n_B = 8
    objSrouby.d_B0 = 16                                                                                                       # jmenovity prumer zavitu sroubu                                [mm]
    objSrouby.d_Bs = objSrouby.d_B0 
    objSrouby.d_B4 = 24
    objSrouby.p_t = 2                                                                                                         # stoupani zavitu                                               [mm]
    objSrouby.l_B = 35                                                                                                        # obr3 - delka zatizene casti sroubu                                   [mm]
    objSrouby.f_B0 = 500                                                                                                      # jmenovite (dovolene) napeti ve sroubu                         [MPa]
    objSrouby.Eps1_plus = 0
    objSrouby.Eps1_minus = 0
    objSrouby.A = 13
    objSrouby.mu_t = 0.2
    objSrouby.mu_n = 0.2
    objSrouby.alpha = 30


    objSrouby.T_Ezk = numpy.asarray([20,100])
    objSrouby.E_zk = numpy.asarray([195000,195000])
    objSrouby.T_azk = numpy.asarray([20,100])
    objSrouby.alfa_zk = numpy.asarray([12e-6,12e-6])

    ##### objTesneni:  ######
    objTesneni = VolbaTesneni(1) 
    objTesneni.druh = 2         # 1-kov, 2-nekov
    # Pro vsechny typy tesneni:
    objTesneni.Q_A = 10       # priloha G - neni pozadovana mira netesnosti   [MPa]
    objTesneni.d_G1 = 110       # teoreticky vnitrni prumer tesnici plochy      [mm]
    objTesneni.d_G2 = 150      # teoreticky vnejsi prumer tesnici plochy       [mm]
    objTesneni.e_G = 2         # tloustka tesneni v nezatizenem stavu          [mm]
    objTesneni.Q_smax = 120    # maximalni dovoleny tlak na tesneni            [MPa]

    objTesneni.mu_G = 0.05
    # Zkouska P_QR
    objTesneni.T_PQR = numpy.asarray([20,23,150])

    objTesneni.Q_I = numpy.asarray([1,30,50])       # pocatecni napeti v tesneni                    [MPa]
    objTesneni.Q_R = numpy.asarray([1,30,50])       # zbytkove naapeti v tesneni                    [MPa]
    objTesneni.d_Gext = numpy.asarray([92,92,92])   # vnejsi prumer tesneni pouziteho pri zkousce   [mm]
    objTesneni.d_Gint = numpy.asarray([49,49,49])   # vnitrni prumer tesneni pouziteho pri zkousce  [mm]
    objTesneni.K = numpy.asarray([500000,500000,500000])    # tuhost zk. zarizeni                           [N/mm]


    # Zkouska tloustky tesneni
    objTesneni.T_Gzk = numpy.asarray([20,160])
    objTesneni.Q_Gzk = [numpy.asarray([5.07,7.57,10.07,12.56,15.05,17.55,20.05,30.02,39.99,49.98,59.96,79.88]),numpy.asarray([5.07,7.57,10.07,12.56,15.05,17.55,20.05,30.02,39.99,49.98,59.96,79.88])]
    objTesneni.e_Gzk = [numpy.asarray([0.0119,0.0204,0.0296,0.0402,0.0519,0.0648,0.0789,0.1735,0.3031,0.4111,0.4958,0.6145]),numpy.asarray([0.0119,0.0204,0.0296,0.0402,0.0519,0.0648,0.0789,0.1735,0.3031,0.4111,0.4958,0.6145])]

    #objTesneni.Q_Gzk = [numpy.asarray([5.07,7.57,10.07,12.56,15.05,17.55,20.05,30.02,39.99,49.98,59.96,79.88]),numpy.asarray([1.93,2.53,3.52,4.40,5.20,5.99,6.88,10.15,13.42,16.80,20.08])]
    #objTesneni.e_Gzk = [numpy.asarray([0.0119,0.0204,0.0296,0.0402,0.0519,0.0648,0.0789,0.1735,0.3031,0.4111,0.4958,0.6145]),numpy.asarray([0.0109,0.0815,0.1781,0.3119,0.4428,0.5525,0.6388,0.7117,0.9539,1.0841,1.1692,1.2312])]

    # Zkouska modulu pruznosti
    objTesneni.T_Ezk = numpy.asarray([20,160])
    objTesneni.Q_Ezk = [numpy.asarray([5.07,7.57,10.07,12.56,15.05,17.55,20.05,30.02,39.99,49.98,59.96,79.88]),numpy.asarray([5.07,7.57,10.07,12.56,15.05,17.55,20.05,30.02,39.99,49.98,59.96,79.88])]
    objTesneni.E_Ezk = [numpy.asarray([1120,1210,1194,1231,1284,1330,1396,1733,2032,2328,2664,3305]),numpy.asarray([1120,1210,1194,1231,1284,1330,1396,1733,2032,2328,2664,3305])]

    #objTesneni.Q_Ezk = [numpy.asarray([5.07,7.57,10.07,12.56,15.05,17.55,20.05,30.02,39.99,49.98,59.96,79.88]),numpy.asarray([1.93,2.53,3.52,4.40,5.20,5.99,6.88,10.15,13.42,16.80,20.08])]
    #objTesneni.E_Ezk = [numpy.asarray([1120,1210,1194,1231,1284,1330,1396,1733,2032,2328,2664,3305]),numpy.asarray([94,90,104,132,176,220,218,213,315,413,390])]

	# Zkouska teplotni roztaznosti
    objTesneni.T_azk = numpy.asarray([20,100])
    objTesneni.alfa_zk = numpy.asarray([7.9e-6,7.9e-6])



    ##### objZatizeni:  ######
    objTesneni.T = numpy.asarray([objTesneni.T0,40,60,80,100,120,140])
    lenP = len(objTesneni.T)
    objTesneni.Q_sminLI = numpy.asarray([None,5,5,5,5,5,5])    # minimalni povrchovy (utahovaci) tlak          [MPa]
                        # pusobici na tesneni , pozadovany pro tridu tesnosti L v podminkach provozu
    objSrouby.T = objTesneni.T
    objPrvniPriruba.T = objTesneni.T
    objDruhaPriruba.T = objPrvniPriruba.T
    objPrvniPodlozka.T = objTesneni.T
    objPrvniPodlozka.alfa = numpy.full((lenP),0)
    objDruhaPodlozka = copy.deepcopy(objPrvniPodlozka)


    objZatizeni.F_ZI = numpy.full((lenP),0)
    objZatizeni.P_I = numpy.full((lenP),0)
    objZatizeni.F_XI = numpy.full((lenP),0)
    objZatizeni.F_YI = numpy.full((lenP),0)
    objZatizeni.M_XI = numpy.full((lenP),0)
    objZatizeni.M_YI = numpy.full((lenP),0)
    objZatizeni.M_ZI = numpy.full((lenP),0)
    objZatizeni.N_R = 1


    ##### objMatice:  ######
    objMatice.e_N = 13
    objMatice.f_N = 800

    ###########################################################################
    ##############################  VYPOCET  ##################################
    ###########################################################################

    objPrvniPriruba.setn_B(objSrouby.n_B)
    objPrvniPriruba.sete()
    objDruhaPriruba.sete()
    objDruhaPriruba.n_B = objPrvniPriruba.n_B
    objTesneni.sete()

    # Prvni dilci vypocty
    objPrvniPriruba.VypocitejPrirubu()
    objDruhaPriruba.VypocitejPrirubu()
    objSrouby.VypocitejSrouby()
    objPrvniPodlozka.calc635(objPrvniPriruba.d_5,objSrouby.d_B4,objPrvniPriruba.n_B)
    objDruhaPodlozka.calc635(objDruhaPriruba.d_5,objSrouby.d_B4,objDruhaPriruba.n_B)

    # Parametry tesneni
    objTesneni.calc642()

    objZatizeni.setTypVypoctu(TypVypoctu)
    objZatizeni.setall(objPrvniPriruba, objDruhaPriruba, objSrouby, objTesneni, objMatice, objPrvniPodlozka, objDruhaPodlozka)
    
    objZatizeniA = copy.deepcopy(objZatizeni)
    objZatizeniB = copy.deepcopy(objZatizeni)

    objZatizeniA.solve(+1)
    objZatizeniB.solve(-1)

    if objZatizeniA.F_B0nom >= objZatizeniB.F_B0nom:
        objZatizeni = copy.deepcopy(objZatizeniA)
    else:
        objZatizeni = copy.deepcopy(objZatizeniB)

    # Pomery zatizeni
    objZatizeni.calc8()

    neniSplnenaPodminka = objZatizeni.conditionl_B() ## dodelat hlasku
    if neniSplnenaPodminka:
        print('Neni splnena podminka delky (98)!')
        #ZDE - prerusi vypocet, vypise hlasku, neshodi cely "program"


    print('Done!')
if __name__ == '__main__':
    sys.exit(int(main() or 0))

    
