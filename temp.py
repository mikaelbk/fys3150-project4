from numpy import *
from matplotlib.pyplot import *
rc('text',usetex=True)

x = linspace(0,100,101)
y = sin(x)

plot(x,y)
show()