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

def main():
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
            10: TocivaPrirubaSObrubou_Lemem(0,2),
            11: ObecnaPriruba(1,1,1),
            12: TocivaPrirubaSObrubou_Lemem(1,0),
            }[typ]

    objPrvniPriruba = VolbaPriruby(7)
    objPrvniPriruba.alfa = numpy.asarray([11.0e-6,11.0e-6])
    objPrvniPriruba.d_0 = 139.8
    objPrvniPriruba.d_3 = 200
    objPrvniPriruba.d_4 = 235
    objPrvniPriruba.d_5 = 19
    objPrvniPriruba.d_S = 142.8
    objPrvniPriruba.E = numpy.asarray([190000,190000])
    objPrvniPriruba.e_1 = 3
    objPrvniPriruba.e_F = 16
    objPrvniPriruba.e_Fb = 16
    objPrvniPriruba.e_Ft = 16
    objPrvniPriruba.sete()

    objDruhaPriruba = VolbaPriruby(7)
    objPrvniPriruba.alfa = numpy.asarray([11.0e-6,11.0e-6])
    objDruhaPriruba.d_0 = objPrvniPriruba.d_0
    objDruhaPriruba.d_3 = objPrvniPriruba.d_3
    objDruhaPriruba.d_4 = objPrvniPriruba.d_4
    objDruhaPriruba.d_5 = objPrvniPriruba.d_5
    objDruhaPriruba.d_S = 143.3
    objDruhaPriruba.E = objPrvniPriruba.E
    objDruhaPriruba.e_1 = 3.5
    objDruhaPriruba.e_F = objPrvniPriruba.e_F
    objDruhaPriruba.e_Fb = objPrvniPriruba.e_Fb
    objDruhaPriruba.e_Ft = objPrvniPriruba.e_Ft
    objDruhaPriruba.sete()

    objSrouby = Sroub(1)
    objSrouby.E = numpy.asarray([190000,190000])
    objSrouby.alfa = numpy.asarray([11.0e-6,11.0e-6])
    objSrouby.d_B0 = 16
    objSrouby.d_B4 = 24
    objSrouby.l_B = 33.5
    objSrouby.n_B = 8
    objSrouby.p_t = 2

    def VolbaTesneni(typ):
        A = {
            1 : TesneniTyp1,
            2 : TesneniTyp2,
            3 : TesneniTyp3,
            4 : TesneniTyp4,
            }[typ]
        return A()

    objTesneni = VolbaTesneni(1)
    objTesneni.d_G1 = 145
    objTesneni.d_G2 = 182
    objTesneni.e_G = 1.5
    objTesneni.sete()

    objMatice = Matice()
    objMatice.e_N = 13

    objPrvniPodlozka = Podlozka(0)
    objDruhaPodlozka = Podlozka(0)
    objPrvniPodlozka.e = 0
    objPrvniPodlozka.E = numpy.asarray([205000,205000])
    objDruhaPodlozka.e = 0
    objDruhaPodlozka.E = numpy.asarray([205000,205000])
    objZatizeni = Zatizeni()
    objPrvniPriruba.setn_B(objSrouby.n_B)
    objDruhaPriruba.n_B = objPrvniPriruba.n_B
    objTesneni.E = numpy.asarray([1284,1284])
    objTesneni.alfa = numpy.asarray([1.2e-4,1.2e-4])

    # Prvni dilci vypocty
    objPrvniPriruba.VypocitejPrirubu()
    objDruhaPriruba.VypocitejPrirubu()
    objSrouby.VypocitejSrouby()
    objPrvniPodlozka.calc635(objPrvniPriruba.d_5,objSrouby.d_B4,objPrvniPriruba.n_B)
    objDruhaPodlozka.calc635(objDruhaPriruba.d_5,objSrouby.d_B4,objDruhaPriruba.n_B)

    # Parametry tesneni
    objTesneni.calc642()

    # Volba typu vypoctu
    def ZvolTypVypoctu(typ):
        return {
            1 : "MiraNetesnosti",
            2 : "KontrolaSroubu",
            }[typ]
    TypVypoctu = ZvolTypVypoctu(1)

    objZatizeni.setTypVypoctu(TypVypoctu)
    objZatizeni.setall(objPrvniPriruba, objDruhaPriruba, objSrouby, objTesneni, objMatice, objPrvniPodlozka, objDruhaPodlozka)
    objZatizeni.setF_G0()

    # Prvni aproximace
    objTesneni.calc643(objPrvniPriruba,objDruhaPriruba,objZatizeni.F_G0)

    # Iterace tesneni
    objTesneni.iteraceb(objPrvniPriruba,objDruhaPriruba,objZatizeni.F_G0)    
      
    objZatizeni.calc7()
    
    # Iterace sily
    objZatizeni.iteraceF()

    # 7.5.2 Zohledneni rozptylu zatizeni sroubu pri montazi
    objZatizeni.calc752()
    
    # Vyhodnoceni
    objZatizeni.calc76()
    objZatizeni.calc8()

    neniSplnenaPodminka = objZatizeni.conditionl_B() ## dodelat hlasku
    if neniSplnenaPodminka:
        print('Neni splnena podminka delky (98)!')
        sys.exit(int(0))

    # Pomery zatizeni
    objZatizeni.calcPhi_B()
    objZatizeni.calcPhi_G()





    objZatizeni.calcF_GI()

    objZatizeni.calcPhi_G()
    Theta_F1 = objZatizeni.calcTheta_F(objPrvniPriruba)
    Theta_F2 = objZatizeni.calcTheta_F(objDruhaPriruba)

    #F_G0 = objSrouby.calcF_B0req(objZatizeni.F_G0req)
    #F_G0req = objTesneni.calcF_G0req(objZatizeni.F_G0)
    #objTesneni.F_G0req
    #print(F_G0)
    #print(F_G0req)


    objPrvniPriruba.calch_QGHL(objTesneni.d_Ge)
    objDruhaPriruba.calch_QGHL(objTesneni.d_Ge)

    objZatizeni.calcPhi_G()
    print('Done!')
if __name__ == '__main__':
    sys.exit(int(main() or 0))

    
