from Soucast import *
from math import *

class Podlozka(Soucast):
    """description of class"""

    d_W1 = 24
    d_W2 = 30
    pouzito = True

    def calcX_W(self,d_5,d_B4,n_B):
        """(43)(44)(45)(46)(47)(48)(50)"""
        if self.pouzito:
            self.b_W = (self.d_W2 - self.d_W1) / 2
            self.d_W = (self.d_W2 + self.d_W1) / 2
            self.d_K1 = max(d_5,self.d_W1)
            self.d_K2 = min(d_B4,self.d_W2)
            self.b_KB = (self.d_K2 - self.d_W1) / 2
            self.X_W = self.e / (n_B * pi * self.d_W * self.b_W) * (2 * self.b_W / (self.b_W + self.b_KB) \
                + self.e / (self.b_W - self.b_KB)) / (1 + self.e / (self.b_W - self.b_KB))
        else:
            X_W = 0