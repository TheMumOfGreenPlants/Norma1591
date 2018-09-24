from math import pi, acos

class Zatizeni(object):
    """description of class"""
    P_I = 0
    F_XI = 0
    F_YI = 0
    F_ZI = 0
    M_XI = 0
    M_YI = 0
    M_ZI = 0

    T_0 = 0
    T_BI = 0
    T_1FI = 0
    T_2FI = 0
    T_1LI = 0
    T_2LI = 0
    T_1WI = 0
    T_2WI = 0
    T_GI = 0

    alfa_BI = 0
    alfa_1FI = 0
    alfa_2FI = 0
    alfa_1LI = 0
    alfa_2LI = 0
    alfa_1WI = 0
    alfa_2WI = 0
    alfa_GI = 0
        
    def calcA_Q(self, d_Gi):
        """(90)"""
        self.A_Q = (pi * d_Gi**2) /4

    def calcF_QI(self, d_Gi):
        """(91)"""
        self.calcA_Q(d_Gi)
        self.F_IQ = self.A_Q * P_I

    def calcFM(self, d_3e, varianta):
        """(92)(93)(94)(95)(96)"""
        self.F_AI = F_ZI
        self.F_LI = (F_XI**2 + F_YI**2)**(1/2)
        self.M_AI = (M_XI**2 + M_YI**2)**(1/2)
        self.M_TGI = M_ZI
        def f(var):
            return {
                a : F_AI + (4 / d_3e) * M_AI,
                b : F_AI - (4 / d_3e) * M_AI,
                }[var]
        self.F_RI = f(varianta)

    def calcdeltaU_TI(self, e_Ft1, e_Ft2, e_L1, e_L2, e_Gt, e_W1, e_W2, l_B):
        """(97)"""
        self.deltaU_TI = l_B * self.alfa_BI * (self.T_BI - self.T_0) \
                         - self.e_Ft1 * self.alfa_1FI * (self.T_1FI - T_0) \
                         - self.e_Ft1 * self.alfa_1FI * (self.T_1FI - T_0) \
                         - self.e_Ft2 * self.alfa_2FI * (self.T_2FI - T_0) \
                         - self.e_Lt1 * self.alfa_1LI * (self.T_1LI - T_0) \
                         - self.e_Lt2 * self.alfa_2LI * (self.T_2LI - T_0) \
                         - self.e_W1 * self.alfa_1WI * (self.T_1WI - T_0) \
                         - self.e_W2 * self.alfa_2WI * (self.T_2WI - T_0) \
                         - self.e_Gt * self.alfa_2WI * (self.T_2WI - T_0)

    def conditionl_B(self, e_Ft1, e_Ft2, e_L1, e_L2, e_Gt, e_W1, e_W2, l_B):
        """(98)"""
        sum = e_Ft1 + e_Ft2 + e_L1 + e_L2 + e_Gt + e_W1 + e_W2
        return sum == l_B

    def calcY(self):
        """(99)(100)(101)(102)"""

        
