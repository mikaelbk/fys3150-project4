from matplotlib.pyplot import *
from metropolis_ising import *

L = 40
cycles = int(1000)
T = 1.75
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

meanE[0] = energy[0]
for i in range(1,len(energy)):
	meanE[i] = sum(energy[:i])/i

print 'plotting'
counter = 0
for i in range(len(energy)):
	if(meanE[i] >= counter * meanE[-1] / 8):
		subplot(3,3,counter+1)
		imshow(s[i],cmap = 'Greys', interpolation = 'nearest')
		percentage = str((100.*i)/len(energy))
		title('MCC# ' + str(i) + '. E = ' + str(energy[i]/L**2))
		counter = counter + 1
		axis('off')
show()
plot(energy/L**2, label = "microstate energy")
plot(meanE/L**2, label = "system average energy")
legend(loc = 'best')
xlabel("Number of montecarlo cycles")
ylabel("Energy in kT/J")
axis([0,1000,0,2.5])
show()
#plot(meanE[int(9*len(meanE)/10.):] - L2expectE(1))