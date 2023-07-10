import numpy as np
import random 
import matplotlib.pyplot as plt
def analytic(t,N,e):
    return 1/2*( np.exp(-t)+np.exp(t)*np.exp(-N*e))


class Lattice(object):
    def __init__(self, N, a, N_thermal, N_mc, N_correlation, Delta) -> None:
        self.N_thermal = N_thermal
        self.N_correlation = N_correlation
        self.N_mc = N_mc
        self.a = a
        self.N = N
        self.Delta = Delta
        self.lattice = np.zeros(N)
        self.total_time = N*a
        self.sweeps = 0
        self.thermalize()
    def sweep(self):
        for j in range(self.N):
            pos =j
            delta = random.uniform(-self.Delta,self.Delta)
            candidate = self.lattice[pos] + delta
            next_pos = (pos+1)%(self.N)
            pre_pos = (pos-1)%(self.N)
            dS = 1/2* ((self.lattice[next_pos]-candidate)**2/self.a) + 1/2 *((candidate-self.lattice[pre_pos])**2/self.a) -1/2 *((self.lattice[next_pos]-self.lattice[pos])**2/self.a) -1/2* ((self.lattice[pos]-self.lattice[pre_pos])**2/self.a) + self.a*1/8 *(self.lattice[next_pos]+candidate)**2 +self.a*1/8 *(candidate+self.lattice[pre_pos])**2-self.a*1/8 *(self.lattice[next_pos]+self.lattice[pos])**2-self.a*1/8 *(self.lattice[pre_pos]+self.lattice[pos])**2
            probs = random.random()
            if probs < min(1,np.exp(-dS)):
                self.lattice[pos] = candidate
    def thermalize(self):
        for i in range(self.N_thermal):
            self.sweep()
    def measure_twopoint(self):
        positions = [i for i in range(self.N)]
        values = [[] for i in range(len(positions))]
        for i in range(self.N_mc):
            print(i)
            self.sweep()
            if self.sweeps % self.N_correlation ==0:
                for i in range(len(positions)):
                    values[i] += [self.lattice[0]*self.lattice[positions[i]]]
            self.sweeps+=1
        means = [np.mean(value) for value in values] + [np.mean(values[0])]
        times = [self.a*position for position in positions] + [self.a*self.N]
        return means,times



def main():
    N = 16
    a = 0.3125
    lat = Lattice(N=N,a=a,N_thermal=100,N_mc=10**5,N_correlation=3,Delta=1)
    means,times = lat.measure_twopoint()
    
    ts = np.linspace(start=0,stop=N*a)
    corr = [analytic(t,N,a) for t in ts]
    plt.plot(times,means,'o')
    plt.plot(ts,corr)
    plt.show()
main()