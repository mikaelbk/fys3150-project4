from numpy import *
from matplotlib.pyplot import *

def expectE(T):
	b = 1/T
	return -8*( e**(8*b) - e**(-8*b) ) / ( e**(8*b) + e**(8*b) + 6 )

print(expectE(1)/4)