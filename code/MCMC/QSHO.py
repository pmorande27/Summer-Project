import numpy as np
import random 
import matplotlib.pyplot as plt

def main():
    # Parameters
    N = 1000
    e = 0.01
    N_sweeps = 100
    Delta = 0.5
    total_time = N*e
    sites = np.zeros(N)
    for i in range(N_sweeps):
        for j in range(N):
            pos = random.randint(0,N-1)
            delta = random.uniform(-Delta,Delta)
            candidate = sites[pos] + delta
            next_pos = (pos+1)%N
            pre_pos = (pos-1)%N
            dS = 1/2* ((sites[next_pos]-candidate)/e)**2 + 1/2 *((candidate-sites[pre_pos])/e)**2 -1/2 *((sites[next_pos]-sites[pos])/e)**2 -1/2* ((sites[pos]-sites[pre_pos])/e)**2 + 1/2 *candidate**2 - 1/2* sites[pos]**2
            probs = random.random()
            if probs < min(1,np.exp(-dS)):
                sites[pos] = candidate
    plt.plot(sites)
    plt.show()
main()