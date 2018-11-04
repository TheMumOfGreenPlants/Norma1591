from PrirubaSKuzelovymKrkem import *
from Podlozka import *
from Sroub import *
from Tesneni import *
from Zatizeni import *
from math import pi, acos
from sys import stdin
import sys

def main():
    objPrvniPriruba = PrirubaSKuzelovymKrkem()
    objPrvniPriruba.E = 200000
    objPrvniPriruba.alfa = 11.3e6
    objDruhaPriruba = PrirubaSKuzelovymKrkem()
    objDruhaPriruba.E = 200000
    objDruhaPriruba.alfa = 11.3e6
    objSrouby = Sroub()
    objTesneni = Tesneni()
    objMatice = Matice()
    objPrvniPodlozka = Podlozka()
    objDruhaPodlozka = Podlozka()
    objZatizeni = Zatizeni()
    objDruhaPriruba.e_1 = 11.1
    objDruhaPriruba.e_S = 11.1
    objDruhaPriruba.l_H = 42
    objDruhaPriruba.d_1 = 70.55
    objPrvniPriruba.getn_B(objSrouby.n_B)
    objDruhaPriruba.getn_B(objSrouby.n_B)
    objTesneni.E = 2103
    objTesneni.alfa = 16.4e6

       
    objSrouby.calcA_B()
    objTesneni.setPriruby(objPrvniPriruba, objDruhaPriruba)
    objPrvniPriruba.calcZ_FL()
    objDruhaPriruba.calcZ_FL()
    objPrvniPriruba.calcd_3e()
    objDruhaPriruba.calcd_3e()
    objSrouby.calcF_R0()

    objZatizeni.conditionl_B(objPrvniPriruba.e_Ft,objDruhaPriruba.e_Ft,objPrvniPriruba.e_L,objDruhaPriruba.e_L,objPrvniPodlozka.e_W,objDruhaPodlozka.e_W,objTesneni.e_G,objSrouby.l_B)

    F_G0 = 282018.6  # F_G0pocatecni - vlastni volba 
    F_G0 = objSrouby.A_B * objSrouby.f_B0 / 3 - objSrouby.F_R0
    F_G0req = 0
    while abs(F_G0req - F_G0) >= (F_G0req * 0.001):
        objTesneni.calcF_G0req(F_G0)
        objTesneni.F_G0req
        objSrouby.calcPreload(objTesneni.F_G0req)
        F_G0 = objSrouby.calcF_B0req(F_G0req)
        F_G0req = objTesneni.calcF_G0req(F_G0)
        objTesneni.F_G0req
        print(F_G0)
        print(F_G0req)


    objPrvniPriruba.calch_QGHL(objTesneni.d_Ge)
    objDruhaPriruba.calch_QGHL(objTesneni.d_Ge)
    objSrouby.calcX_B()
    vysledek = objSrouby.Preload

    objZatizeni.setall(objPrvniPriruba, objDruhaPriruba, objSrouby, objTesneni, objMatice)
    objZatizeni.calcY()
    objTesneni.calcA_Gt()
    objZatizeni.calcPhi_G()
    print('Done!')
if __name__ == '__main__':
    sys.exit(int(main() or 0))

    
