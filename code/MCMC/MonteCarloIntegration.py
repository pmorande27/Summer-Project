import random
import numpy as np
import matplotlib.pyplot as plt
def integrable(x):
    return x*(1-x)
def montecarloint(N):
    b = 1
    a = 0
    results =[0 for i in range(N)]
    for i in range(N):
        sample = random.random()
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
    Nset = [10**3,5*10**3,10**4,10**5,10**6,10**7,10**8]
    results = [montecarloint(N) for N in Nset ]
    ys = [0]*len(Nset)
    es = [0]*len(Nset)
    for i in range(len(results)):
        ys[i],es[i] = get_value(results[i])
    plt.errorbar(Nset, ys, yerr=es, fmt='s')
    plt.axhline(y = 0.16666666666666666, color = 'r',linestyle='dotted')
    print(ys)
    plt.xscale('log')
    plt.ylim(0.15,0.18)
    plt.show()

plot()
