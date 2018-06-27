from PrirubaSKuzelovymKrkem import *
from Sroub import *
from Tesneni import *
from Zatizeni import *
from math import pi, acos
from sys import stdin
import sys

def main():
    objPrvniPriruba = PrirubaSKuzelovymKrkem()
    objDruhaPriruba = PrirubaSKuzelovymKrkem()
    objSrouby = Sroub()
    objTesneni = Tesneni()
    objZatizeni = Zatizeni()
    objDruhaPriruba.e_1 = 11.1
    objDruhaPriruba.l_H = 42
    objDruhaPriruba.d_1 = 70.55
       
    objSrouby.calcA_B(objPrvniPriruba.n_B)
    objTesneni.setPriruby(objPrvniPriruba, objDruhaPriruba)
    objPrvniPriruba.calcZ_FL()
    objDruhaPriruba.calcZ_FL()
    objPrvniPriruba.calcd_3e()
    objDruhaPriruba.calcd_3e()
    objSrouby.calcF_R0()

    F_G0 = 282018.6  # F_G0pocatecni - vlastni volba 
    F_G0 = objSrouby.A_B * objSrouby.f_b0 / 3 - objSrouby.F_R0
    F_G0req = 773632.2
    while (F_G0req - F_G0) >= (F_G0req * 0.001):
        objTesneni.calcF_G0req(F_G0)
        objTesneni.F_G0req
        objSrouby.calcPreload(objTesneni.F_G0req)
        F_G0 = objSrouby.calcF_B0req(F_G0req)
        F_G0req = objTesneni.calcF_G0req(F_G0)
        objTesneni.F_G0req
        print(F_G0)
        print(F_G0req)

    vysledek = objSrouby.Preload
if __name__ == '__main__':
    sys.exit(int(main() or 0))

    
