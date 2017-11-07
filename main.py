from matplotlib.pyplot import *
from metropolis_ising import metropolis_ising

N = 10
for i in range(N):
	if(i<N/2):
		s = metropolis_ising()
		filename = 'E0_' + str(i) + '_.png'
	else:
		s = metropolis_ising(rand = True)
		filename = 'random_' + str(i) + '_.png'
	imshow(s,cmap = 'Greys', interpolation = 'nearest')
	savefig(filename)
