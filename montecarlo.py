from matplotlib.pyplot import *
from metropolis_ising import *

L = 30
cycles = int(5E2)
s = zeros((cycles,L,L))
s[0] = metropolis_init(rand = True)
energy = zeros(cycles)
energy[0] = totalEnergy(s[0])

for i in range(cycles-1):
	s[i+1] = metropolis_ising(s[i],T=1.5)
	energy[i+1] = totalEnergy(s[i+1])

for i in range(4):
	subplot(2,2,int(i+1))
	imshow(s[int(float(i)*(len(s)-1)/4.)],cmap = 'Greys', interpolation = 'nearest')
show()
plot(energy)
show()

#plot(energy)
#show()