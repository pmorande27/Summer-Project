import numpy as np
import random
import matplotlib.pyplot as plt
def integrable(x):
    return (1-np.exp(-3))/(1+x/9)
def inverse(sampleU):
    return -np.log(1-(1-np.exp(-3))*sampleU)
def montecarloint(N):
    b = 3
    a = 0
    results =[0 for i in range(N)]
    for i in range(N):
        sampleU = random.uniform(0,1)
        sample = inverse(sampleU)
        results[i] = integrable(sample)
    return results

def get_value(results):
    values = np.array(results)
    result = np.mean(values)
    sq_values  = np.square(values)
    average_sq_values = np.mean(sq_values)
    uncertainty = np.sqrt((average_sq_values-result**2)/len(values))
    return(result,uncertainty)

def plot():
    Nset = [10**2,10**3,5*10**3,10**4,10**5]
    results = [montecarloint(N) for N in Nset ]
    ys = [0]*len(Nset)
    es = [0]*len(Nset)
    for i in range(len(results)):
        ys[i],es[i] = get_value(results[i])
    plt.errorbar(Nset, ys, yerr=es, fmt='s')
    plt.axhline(y = 0.873109, color = 'r',linestyle='dotted')
    print(ys)
    plt.xscale('log')
    plt.ylim(0.8,1)
    plt.show()

plot()

