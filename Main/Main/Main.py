from PrirubaSKuzelovymKrkem import *
from IntegralniPriruba import *
from Podlozka import *
from Sroub import *
from Tesneni import *
from Zatizeni import *
from math import pi, acos
from sys import stdin
import sys

def main():
    objPrvniPriruba = PrirubaSKuzelovymKrkem()
    objPrvniPriruba.E = numpy.asarray([200000,200000])
    objPrvniPriruba.alfa = numpy.asarray([11.3e-6,11.3e-6])
    objPrvniPriruba.sete()
    objDruhaPriruba = PrirubaSKuzelovymKrkem()
    objDruhaPriruba.E = numpy.asarray([200000,200000])
    objDruhaPriruba.alfa = numpy.asarray([11.3e-6,11.3e-6])
    objDruhaPriruba.sete()
    objSrouby = Sroub()
    objSrouby.E = numpy.asarray([205000,205000])
    objSrouby.alfa = numpy.asarray([11.8e-6,11.8e-6])
    objTesneni = Tesneni()
    objTesneni.sete()
    objMatice = Matice()
    objPrvniPodlozka = Podlozka()
    objDruhaPodlozka = Podlozka()
    objPrvniPodlozka.e = 0
    objPrvniPodlozka.E = numpy.asarray([205000,205000])
    objDruhaPodlozka.e = 0
    objDruhaPodlozka.E = numpy.asarray([205000,205000])
    objZatizeni = Zatizeni()
    objDruhaPriruba.e_1 = 11.1
    objDruhaPriruba.e_S = 11.1
    objDruhaPriruba.l_H = 42
    objDruhaPriruba.d_1 = 70.55
    objPrvniPriruba.getn_B(objSrouby.n_B)
    objDruhaPriruba.getn_B(objSrouby.n_B)
    objTesneni.E = numpy.asarray([2103,2103])
    objTesneni.alfa = numpy.asarray([16.4e-6,16.4e-6])

    obj = IntegralniPriruba()
    obj.calcc_M()

    
    objSrouby.calcA_B()
    objTesneni.setPriruby(objPrvniPriruba, objDruhaPriruba)
    objPrvniPriruba.calcZ_FL()
    objDruhaPriruba.calcZ_FL()
    objPrvniPriruba.calcd_3e()
    objDruhaPriruba.calcd_3e()

    objZatizeni.setall(objPrvniPriruba, objDruhaPriruba, objSrouby, objTesneni, objMatice, objPrvniPodlozka, objDruhaPodlozka)
    neniSplnenaPodminka = objZatizeni.conditionl_B() ## dodelat hlasku
    if neniSplnenaPodminka:
        print('Neni splnena podminka delky (98)!')
        sys.exit(int(0))


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

    
