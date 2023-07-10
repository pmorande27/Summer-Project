import numpy as np
import random 
import matplotlib.pyplot as plt
def analytic(t,N,e):
    return 1/2*( np.exp(-t)+np.exp(t)*np.exp(-N*e))
def main():
    # Parameters
    N = 32
    e = 0.3125
    Nmc = 10**7
    Nsweeps =5
    positions = [i for i in range(32)]
    Delta = 1
    total_time = N*e
    values = [[] for i in range(len(positions))]
    sites = np.zeros(N)
    means = []
    counter =0
    for i in range(Nmc):
        print(i)
        for j in range(N):
            pos =j
            delta = random.uniform(-Delta,Delta)
            candidate = sites[pos] + delta
            next_pos = (pos+1)%N
            pre_pos = (pos-1)%N
            dS = 1/2* ((sites[next_pos]-candidate)**2/e) + 1/2 *((candidate-sites[pre_pos])**2/e) -1/2 *((sites[next_pos]-sites[pos])**2/e) -1/2* ((sites[pos]-sites[pre_pos])**2/e) + e*1/2 *candidate**2 - e*1/2* sites[pos]**2
            probs = random.random()
            if probs < min(1,np.exp(-dS)):
                sites[pos] = candidate
        if counter == Nsweeps:
            for i in range(len(positions)):
                values[i] += [sites[0]*sites[positions[i]]]
            counter =0
        counter+=1
    means = [np.mean(value) for value in values]
    times = [e*position for position in positions]
    ts = np.linspace(start=0,stop=total_time)
    corr = [analytic(t,N,e) for t in ts]
    plt.plot(times,means,'o')
    plt.plot(ts,corr)
    plt.show()
main()