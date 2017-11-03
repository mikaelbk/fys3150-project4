from numpy import *

def metropolis_ising(size = 30., T = 2., rand = False, J=1, b=1.):
	k = 1/(b*T)
	if(rand == 1):
		s = random.choice([-1,1],[size,size])
	else:		
		s = ones([size,size])

	for i in range(100*size**2):
		i = random.randint(size)
		j = random.randint(size)
		top = s[i-1,j]
		left = s[i,j-1]
		bottom = s[i+1,j] if (i != (size-1)) else s[0,j]
		right = s[i,j+1] if (j != (size-1)) else s[i,0]
		Ediff = 2*s[i,j]*(top+bottom+left+right)
		if(Ediff <= 0 or random.random() <= exp(-Ediff/T)):
			s[i,j] *= -1
	return s

def L2expectE():
	return 8*(exp(16*b)-1)/(6*exp(8*b)+exp(16*b)+1)
def L2expectM():
	return 42
def L2expectC():
	return (384*x**2 (e**(-8*b*x) + e**(8*b*x)) + 4)/(k*T**2*(e**(-8*b*x) + e**(8*b*x) + 6)**2)
def L2expectX():
	return 42

if __name__ == '__main__': #only runs if this file is main
	main()

"""
from matplotlib.pyplot import *
imshow(s,cmap = 'Greys', interpolation = 'nearest')
savefig('state_'+ str(k)+'_.png')
"""