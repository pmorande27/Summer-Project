import numpy as np
import random 
import matplotlib.pyplot as plt
from Stats import Stats
def analytic(t, N, e):
    return 1/2*( np.exp(-t)+np.exp(t)*np.exp(-N))


class Lattice(object):
    def __init__(self, N, a, N_thermal, N_measurement, N_sweeps, Delta) -> None:
        self.N_thermal = N_thermal
        self.N_sweeps = N_sweeps
        self.N_measurement = N_measurement
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

            #Entropy change
            dS = (1/(2*self.a)*( (self.lattice[next_pos]-candidate)**2 + (candidate-self.lattice[pre_pos])**2 -(self.lattice[next_pos]-self.lattice[pos])**2 -(self.lattice[pos]-self.lattice[pre_pos])**2 )+ 
            self.a*1/8* ( (self.lattice[next_pos]+candidate)**2 + (candidate+self.lattice[pre_pos])**2-(self.lattice[next_pos]+self.lattice[pos])**2- (self.lattice[pre_pos]+self.lattice[pos])**2))
           
           #Metropolis step
            probs = random.random()
            if probs < min(1, np.exp(-dS)):
                self.lattice[pos] = candidate

    def thermalize(self):

        for i in range(self.N_thermal):
            self.sweep()
    
    def Get_configurations(self):
    
        results = [0 for i in range(self.N_measurement)]
    
        for i in range(self.N_measurement):
    
            for j in range(self.N_sweeps):
    
                self.sweep()
        
            results[i] = self.lattice
    
        return(results)
    
    def generate_measurements(self, observable):

        results = [0 for i in range(self.N_measurement)]

        for i in range(self.N_measurement):

            for j in range(self.N_sweeps):

                self.sweep()

            print(i)
            results[i] = observable(self.lattice)

        return results

    @staticmethod
    def measure_twopoint(lattice):
        N = len(lattice)

        positions = [i for i in range(N)]
        values = [0 for i in range(positions)]
      
        for i in range(len(positions)):
            values[i] = lattice[0]*lattice[positions[i]]

        return values
    
    @staticmethod    
    def measure_greens(lattice):

        N = len(lattice)
        values = [0 for l in range(N)]

        
        for n in range(N):
            values[n] = 1/N* sum([lattice[((j+n)%N)]* lattice[j] for j in range(N)])

        

        return values

    @staticmethod
    def measure_position(lattice):
        return np.average(lattice)
    @staticmethod
    def measure_sq_position(lattice):
        return np.average([x**2 for x in lattice])

def main():
    N = 100
    a = 1
    N_measurements = 10**4
    N_sweeps = 1
    N_thermal = N_sweeps*10
    Delta = 1
    lat = Lattice(N = N, a = a, N_thermal = N_thermal, N_measurement = N_measurements, N_sweeps = N_sweeps, Delta = Delta)
    a = lat.generate_measurements(Lattice.measure_sq_position)
    print(np.average(a))

    

    """ means,times = lat.measure_twopoint()
    
    ts = np.linspace(start = 0, stop = N)
    corr = [analytic(t, N, a) for t in ts]
    plt.plot(times, means, 'o')
    plt.plot(ts, corr)
    #DE = [abs(np.log(means[i]/means[i+1])/a) for i in range(0,N-1)]
    #plt.plot([j for j in range(len(DE))],DE,'o')
    #plt.ylim(0,2)
    plt.show()

    #values, times = lat.measure_acceptance()
    #print(np.mean(values))
    #plt.scatter(times, values)
    #plt.show()
"""
    
main()

