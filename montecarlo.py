from matplotlib.pyplot import *
from metropolis_ising import *


# experimental variables
L = 2
cycles = int(15)
T = 1.0
k = 1.0

# arrays
energy = zeros(cycles)
meanE = zeros(cycles)
s = zeros((cycles,L,L))

# initializing
s[0] = metropolis_init(L = L,rand = True)
energy[0] = totalEnergy(s[0])

# sampling
print 'MC with ' + str(cycles) + ' cycles'
for i in range(cycles-1):
	s[i+1] = metropolis_ising(s[i],T=T,k=k)
	energy[i+1] = totalEnergy(s[i+1])

print 'mean energy'
# mean energy
meanE[0] = energy[0]
for i in range(1,len(energy)):
	meanE[i] = sum(energy[:i])/i
print 'check'

# plotting
"""
counter = 0
for i in range(len(energy)):
	if(meanE[i] >= counter * meanE[-1] / 8):
		subplot(3,3,counter+1)
		imshow(s[i],cmap = 'Greys', interpolation = 'nearest')
		title('MC cycle #' + str(i) + ' (' + str((100.*i)/len(energy)) + '  %)')
		counter = counter + 1
show()
"""

plot(meanE[int(9*len(meanE)/10.):] - L2expectE(1))
show()

#plot(energy)
#show()