from numpy import *
from matplotlib.pyplot import *

J = 1

def x(b):
	return e**(8*J*b)

def y(b):
	return e**(-8*J*b)

def expectE(b):
	#return -8*( e**(8*b) - e**(-8*b) ) / ( e**(8*b) + e**(8*b) + 6 )
	return -8*J*(x(b)-y(b)) / (x(b)+y(b)+6)

def expectESquared(b):
	return 64*J*J*(x(b)+y(b)) / (x(b)+y(b)+6)

def expectMabs(b):
	return (4*x(b)+8) / (x(b)+y(b)+6)

print(expectE(1)/4)
print(expectMabs(1)/4)

# <E>
# <|M|>
# <Cv>
# <x>