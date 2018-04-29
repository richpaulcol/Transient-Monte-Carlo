import numpy as np

D = 0.05
C = 100

def HW_DW(D,C):
	return D*(3.320 - 0.021 * C * D**0.01)**(2.173) * np.exp(-0.04125*C*D**0.01)
	
def DW_HW(D,ks):
	return 'not-implemented'
