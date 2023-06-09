import math
import random
import matplotlib.pyplot as plt
import seaborn as sns

from scipy import rand

def s(x):
    return x**2/2
def metropolis(x):
    deltaX = random.uniform(-0.5,0.5)
    xprime = x + deltaX
    deltaS = s(xprime)-s(x)
    r = random.random()
    if r < math.exp(-deltaS):
        return xprime
    return x
l = []
#for j in range(100):
x0= 0 
for i in range(10**7):
    x0 = metropolis(x0)
    l.append(x0) 
plt.hist(l, bins = int(18000/15)
             )
plt.show()
