from numpy import *
from matplotlib.pyplot import *
from readtxt import *

t1a = get_data("20x20_t1.000000_aligned")
t1r = get_data("20x20_t1.000000_random")
t2a = get_data("20x20_t2.700000_aligned")
t2r = get_data("20x20_t2.700000_random")

subplot(2,2,1)
plot((t1a["mcc"]),t1a["E"])
subplot(2,2,2)
plot(((t1r["mcc"])),(t1r["E"]))
subplot(2,2,3)
plot(((t2a["mcc"])),(t2a["E"]))
subplot(2,2,4)
plot(((t2r["mcc"])),(t2r["E"]))
show()










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

#print(expectE(1)/4)
#print(expectMabs(1)/4)

# <E>
# <|M|>
# <Cv>
# <x>