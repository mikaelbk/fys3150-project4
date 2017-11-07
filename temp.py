from metropolis_ising import *
from matplotlib.pyplot import *

s = metropolis_init()
s2 = metropolis_ising(s)

imshow(s,cmap = 'Greys', interpolation = 'nearest')
show()
imshow(s2,cmap = 'Greys', interpolation = 'nearest')
show()