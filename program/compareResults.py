from numpy import *
from matplotlib.pyplot import *
from readtxt import *

t1a = get_data("20x20_t1.000000_aligned")
t1r = get_data("20x20_t1.000000_random")
t2a = get_data("20x20_t2.400000_aligned")
t2r = get_data("20x20_t2.400000_random")
"""
figure(1)
subplots_adjust(bottom=-0.2)
subplot(4,1,1)
plot((log10(t1a["mcc"])),t1a["E"])
title('T=1, ordered initial lattice')
axis([0,5,-2.0,-1.996])
#subplot(2,2,2)
#plot(log10(((t1r["mcc"]))),(t1r["E"]))
subplot(4,1,2)
plot(log10(((t2a["mcc"]))),(t2a["E"]))
title('T=2.4, ordered initial lattice')
subplot(4,1,3)
plot(log10(((t1r["mcc"]))),(t1r["E"]), label = 'T=1.0')
plot(log10(((t2r["mcc"]))),(t2r["E"]), label = 'T=2.4')
title('Disordered (random) initial lattice')
legend(loc = 'best')
savefig('4c.pdf')"""

figure(2)
subplot(2,2,1)
t1a_crop = (t1a["E"])[50000::]
hist(t1a_crop,100,histtype="bar", normed=1, color='green')
subplot(2,2,2)
hist((t1r["E"])[50000::], 100, normed = 1)
subplot(2,2,3)
hist((t2a["E"])[50000::], 100, normed = 1)
subplot(2,2,4)
hist((t2r["E"])[50000::], 100, normed = 1)
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