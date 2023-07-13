import numpy as np
import random 
import matplotlib.pyplot as plt
def analytic(t, N, e):
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
        self.accepted = 0
        self.thermalize()
    def sweep(self):
        for j in range(self.N):
            pos =j
            delta = random.uniform(-self.Delta,self.Delta)
            candidate = self.lattice[pos] + delta
            next_pos = (pos+1)%(self.N)
            pre_pos = (pos-1)%(self.N)

            #Entropy change
            dS = (1/(2*self.a)*( (self.lattice[next_pos]-candidate)**2 + (candidate-self.lattice[pre_pos])**2 -(self.lattice[next_pos]-self.lattice[pos])**2 -(self.lattice[pos]-self.lattice[pre_pos])**2 )+ 
            self.a*1/8* ( (self.lattice[next_pos]+candidate)**2 + (candidate+self.lattice[pre_pos])**2-(self.lattice[next_pos]+self.lattice[pos])**2- (self.lattice[pre_pos]+self.lattice[pos])**2))
           
           #Metropolis step
            probs = random.random()
            if probs < min(1, np.exp(-dS)):
                self.lattice[pos] = candidate
                #self.accepted += 1

    def thermalize(self):

        for i in range(self.N_thermal):
            self.sweep()
        self.accepted = 0

    def measure_twopoint(self):

        positions = [i for i in range(self.N)]
        values = [[] for i in range(len(positions))]

        for i in range(self.N_mc):
            self.sweep()
            if self.sweeps % self.N_correlation == 0:
                for i in range(len(positions)):
                    values[i] += [self.lattice[0]*self.lattice[positions[i]]]
            self.sweeps+=1

        means = [np.mean(value) for value in values] + [np.mean(values[0])]
        times = [self.a*position for position in positions] + [self.a*self.N]

        return means,times
    
    def measure_greens(self):

        positions = [i for i in range(self.N)]
        values = [[] for n in range(self.N)]

        for i in range(self.N_mc):
            self.sweep()
            print(i)
            if self.sweeps % self.N_correlation == 0:
                for n in range(self.N):
                    values[n] += [1/self.N* sum([self.lattice[((j+n)%self.N)]* self.lattice[j] for j in range(self.N)])]
            self.sweeps+=1

        means = [np.mean(value) for value in values] + [np.mean(values[0])]
        times = [self.a*n for n in range(self.N)] + [self.a*self.N]

        return means,times

    def measure_acceptance(self):
        values = np.zeros(self.N_mc)
        for i in range(self.N_mc):
            self.sweep()
            self.sweeps+=1
            values[i] = float(self.accepted)/self.N
            self.accepted = 0
        return values,[i for i in range(1,self.N_mc+1)]



def main():
    N = 100
    a = 0.05
    N_mc = 10**6
    N_correlation = 20
    N_thermal = 200
    Delta = 1.4
    lat = Lattice(N = N, a = a, N_thermal = N_thermal, N_mc = N_mc, N_correlation = N_correlation, Delta = Delta)

    means,times = lat.measure_greens()
    
    ts = np.linspace(start = 0, stop = N*a)
    corr = [analytic(t, N, a) for t in ts]
    #plt.plot(times, means, 'o')
    #plt.plot(ts, corr)
    DE = [abs(np.log(means[i]/means[i+1])/a) for i in range(0,N-1)]
    plt.plot([j for j in range(len(DE))],DE,'o')
    plt.ylim(0,2)
    plt.show()

    #values, times = lat.measure_acceptance()
    #print(np.mean(values))
    #plt.scatter(times, values)
    #plt.show()

    
main()

