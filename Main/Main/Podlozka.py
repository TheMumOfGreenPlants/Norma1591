from Soucast import *
from math import *

class Podlozka(Soucast):
    """description of class"""

    d_W1 = 0
    d_W2 = 0
    pouzito = True          #mozna nepotrebne

    def calcX_W(self,d_5,d_B4,n_B):
        """(43)(44)(45)(46)(47)(48)(50)"""
        if (self.d_W1 != 0 or self.d_W2 != 0):
            self.b_W = (self.d_W2 - self.d_W1) / 2
            self.d_W = (self.d_W2 + self.d_W1) / 2
            self.d_K1 = max(d_5,self.d_W1)
            self.d_K2 = min(d_B4,self.d_W2)
            self.b_KB = (self.d_K2 - self.d_W1) / 2
            self.X_W = self.e / (n_B * pi * self.d_W * self.b_W) * (2 * self.b_W / (self.b_W + self.b_KB) \
                + self.e / (self.b_W - self.b_KB)) / (1 + self.e / (self.b_W - self.b_KB))
            self.d_n = self.d_K1 + self.b_KB/2
        else:
            self.X_W = 0
            self.d_n = (d_B4 + d_5) / 2