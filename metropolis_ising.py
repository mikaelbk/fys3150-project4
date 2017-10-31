from numpy import *
from matplotlib.pyplot import *

size = 100
T = 2.5

s = random.choice([-1,1],[size,size])

for i in range(100*size**2):
	i = random.randint(size)
	j = random.randint(size)
	top = s[i-1,j]
	left = s[i,j-1]
	bottom = s[i+1,j] if (i != (size-1)) else s[0,j]
	right = s[i,j+1] if (j != (size-1)) else s[i,0]
	Ediff = 2*s[i,j]*(top+bottom+left+right)
	if(Ediff <= 0 or random.random() < exp(-Ediff/T)):
		s[i,j] *= -1

imshow(s,cmap = 'Greys', interpolation = 'nearest')
show()