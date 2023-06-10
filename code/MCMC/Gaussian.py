import math
import random
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import rand
import matplotlib.animation as animation

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
def run(n):
    outputs = np.zeros(n)
    #for j in range(100):
    x0= 0 
    for i in range(n):
        x0 = metropolis(x0)
        outputs[i] = x0 
    return outputs
    """plt.hist(outputs, bins = int(18000/15)
                )
    plt.show()"""
def plot(n):
    plt.hist(run(n), bins = int(18000/15)
                )
    plt.show()
def average(n):
    output = run(n)
    return(np.average(output))

def run_obsevable(n,f):
    outputs = run(n)
    observable = np.zeros(n)
    x0= 0 
    observable[0] = f(x0)
    #for j in range(100):
    for i in range(1,n):
        aux = observable[i-1]
        #print(i)
        observable[i] = (aux*i +f(outputs[i]))/(i+1)
    plt.plot(observable)
    plt.xscale('log')
    plt.ylim(-2,2)
    plt.xlim(1,n)
    plt.show()

run_obsevable(10**6, lambda x: x**3)