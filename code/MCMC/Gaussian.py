import math
import random
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import rand
import matplotlib.animation as animation

def s(x):
    return -np.log(np.exp(-x**2/2)+ np.exp(-(x-20)**2/2))
def metropolis(x,n):
    if n %2 ==0:
        deltaX = random.uniform(-0.5,0.5)
    else: 
        deltaX = random.uniform(-100,100)

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
        x0 = metropolis(x0,i)
        outputs[i] = x0 
    return outputs
    """plt.hist(outputs, bins = int(18000/15)
                )
    plt.show()"""
def plot(n):
    sns.distplot(run(n), hist=False, kde=True, 
             bins=int(180/5), color = 'darkblue', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 4})
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

def error(n,f,number_of_bins):
    outputs = run(n)
    slices = np.split(np.array(outputs),number_of_bins)
    w = len(slices[0])
    aux = sum([ calc_f_w_k(slices,f,number_of_bins,w,k) for k in range(number_of_bins)])
    average_f = aux/number_of_bins
    b = [(sum([f(x) for x in slices[k]])/w - average_f)**2 for k in range(number_of_bins)]
    aux_3 =sum(b)
    aux_3 = np.sqrt(1/(number_of_bins*(number_of_bins-1))*aux_3)
    return aux_3
def calc_f_w_k(slices,f,number_of_bins,w,k):
    aux = 0
    for i in range(number_of_bins):
        if i != k:
            aux += sum([f(x) for x in slices[i]])
    return aux/((number_of_bins-1)*w)

#print(error( 50000, lambda x: x**2, 2500))
plot(10**5)