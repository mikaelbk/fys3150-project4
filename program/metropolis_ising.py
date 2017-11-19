from numpy import *

def energyArray(b = 1):
	array = zeros(17)
	for i in linspace(-8,8,5):
		array[int(i+8)] = exp(-i*b)
	return array

def metropolis_init(L = 20, rand=False):
	if (rand):
		return random.choice([-1,1],[L,L])
	else:
		return ones([L,L])

def metropolis_ising(s, T = 2.0, k = 1.0):
	b = 1/(k*T)
	L = len(s)
	s_out = zeros((L,L))
	energies = energyArray(b)
	for i in range(L):
		for j in range(L):
			s_out[i,j] = s[i,j]
	for a in range(L**2):
		i = random.randint(L)
		j = random.randint(L)
		top = s[i-1,j]
		left = s[i,j-1]
		bottom = s[i+1,j] if (i != (L-1)) else s[0,j]
		right = s[i,j+1] if (j != (L-1)) else s[i,0]
		Ediff = 2*s[i,j]*(top+bottom+left+right)
		if(Ediff <= 0):
			s_out[i,j] *= -1
		elif(random.random() <= energies[int(Ediff+8)]):
			s_out[i,j] *= -1
	return s_out

def totalEnergy(s):
	E = 0
	for i in range(s.shape[0]):
		for j in range(s.shape[1]):
			E = E + s[i,j]*(s[i-1,j]+s[i,j-1])
	return E


def L2expectE(b):
	return 8*(exp(16*b)-1)/(6*exp(8*b)+exp(16*b)+1)
"""
def L2expectM():
	return 42
def L2expectC():
	return (384*x**2 (e**(-8*b*x) + e**(8*b*x)) + 4)/(k*T**2*(e**(-8*b*x) + e**(8*b*x) + 6)**2)
def L2expectX():
	return 42

if __name__ == '__main__': #only runs if this file is main
	main()
"""

"""
from matplotlib.pyplot import *
imshow(s,cmap = 'Greys', interpolation = 'nearest')
savefig('state_'+ str(k)+'_.png')
"""
