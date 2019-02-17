from Priruba import *
from IntegralniPriruba import *
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
            1 : IntegralniPriruba(1,1,1),
            2 : IntegralniPriruba(1,1,2),
            3 : IntegralniPriruba(0,2,1),
            4 : IntegralniPriruba(0,2,2),
            5 : IntegralniPriruba(0,3,1),
            6 : IntegralniPriruba(0,3,2),
            7 : IntegralniPriruba(0,1,1),
            8 : ZaslepovaciPriruba(),
            9 : TocivaPrirubaSObrubou_Lemem(0),
            10: TocivaPrirubaSObrubou_Lemem(0),
            11: IntegralniPriruba(1,1,1),
            12: TocivaPrirubaSObrubou_Lemem(1),
            }[typ]

    objPrvniPriruba = VolbaPriruby(1)
    objPrvniPriruba.E = numpy.asarray([200000,200000])
    objPrvniPriruba.alfa = numpy.asarray([11.3e-6,11.3e-6])
    objPrvniPriruba.sete()

    objDruhaPriruba = VolbaPriruby(1)
    objDruhaPriruba.d_3 = objPrvniPriruba.d_3
    objDruhaPriruba.E = numpy.asarray([200000,200000])
    objDruhaPriruba.alfa = numpy.asarray([11.3e-6,11.3e-6])
    objDruhaPriruba.sete()

    objSrouby = Sroub()
    objSrouby.E = numpy.asarray([210000,210000])
    objSrouby.alfa = numpy.asarray([11.8e-6,11.8e-6])

    def VolbaTesneni(typ):
        A = {
            1 : TesneniTyp1,
            2 : TesneniTyp2,
            3 : TesneniTyp3,
            4 : TesneniTyp4,
            }[typ]
        return A()

    objTesneni = VolbaTesneni(1)
    objTesneni.sete()
    objMatice = Matice()

    objPrvniPodlozka = Podlozka(0)
    objDruhaPodlozka = Podlozka(0)
    objPrvniPodlozka.e = 1
    objPrvniPodlozka.E = numpy.asarray([205000,205000])
    objDruhaPodlozka.e = 1
    objDruhaPodlozka.E = numpy.asarray([205000,205000])
    objZatizeni = Zatizeni()
    #objDruhaPriruba.e_1 = 11.1
    #objDruhaPriruba.e_S = 11.1
    #objDruhaPriruba.l_H = 42
    #objDruhaPriruba.d_1 = 70.55
    objPrvniPriruba.setn_B(objSrouby.n_B)
    objDruhaPriruba.n_B = objPrvniPriruba.n_B
    objTesneni.E = numpy.asarray([1280,1280])
    objTesneni.alfa = numpy.asarray([16.4e-6,16.4e-6])

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

    
