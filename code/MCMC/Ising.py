import numpy as np
import random 
import math
import matplotlib.pyplot as plt
from scipy.ndimage import convolve,generate_binary_structure

class Ising(object):
    def __init__(self, N, kT, N_thermal, N_measuremnt, N_sweeps) -> None:
        self.N = N
        self.N_thermal = N_thermal
        self.N_measurement = N_measuremnt
        self.N_sweeps = N_sweeps
        self.lattice = np.zeros((N,N),dtype=float)
        random = np.random.random((self.N,self.N))
        self.lattice[random>=0.50] = 1
        self.lattice[random<=0.50] = -1
        self.kT = kT
        self.beta = 1/kT
        self.thermalize()


    def __str__(self) -> str:
        """
        Override to string method for printing
        """
        return str(self.lattice)
    
    def sweep(self):
        for x in range(self.N):
            for y in range(self.N):
                i,j = np.random.randint(0,self.N), np.random.randint(0,self.N)
                spin_i = self.lattice[i,j] #initial spin
                spin_f = spin_i*-1 #proposed spin flip
            # compute change in energy
                E_i = 0
                E_f = 0
                E_i += -spin_i*(self.lattice[(i-1) % self.N ,j] + self.lattice[(i+1)%self.N,j] + self.lattice[i,(j-1)%self.N] + self.lattice[i,(j+1) %self.N])
                E_f +=  -spin_f*(self.lattice[(i-1) % self.N ,j] + self.lattice[(i+1)%self.N,j] + self.lattice[i,(j-1)%self.N] + self.lattice[i,(j+1) %self.N])
                delta_E = E_f-E_i
                if delta_E < 0:
                    self.lattice[i,j] = spin_f
                    
                else:
                    if random.random() <= math.e**(-self.beta*delta_E):
                        self.lattice[i,j] = spin_f
    
    @staticmethod
    def get_magnetisation(lattice) -> float:
        """
        Method used to get the magnetisation of the lattice
        """
        return lattice.sum()
    
    @staticmethod
    def get_energy(lattice) -> float:
        """
        Method used to get the energy of the system
        """
        kernel = generate_binary_structure(2,1)
        kernel[1][1] = False
        energy_array = - lattice*convolve(lattice,kernel,mode='constant',cval=0)
        return energy_array.sum()/2
    
    def thermalize(self):

        for i in range(self.N_thermal):
            print(i)
            self.sweep()
    
    def measurement(self, observable):
        results = [0 for i in range(self.N_measurement)]
        for i in range(self.N_measurement):
            for j in range(self.N_sweeps):
                self.sweep()
            print(i)
            results[i] = observable(self.lattice)
        return(results)
    
    def autocorrelation(self, observable):
        steps = [i for i in range(11)]
        measurements = [0 for i in range(self.N_measurement)]
        for i in range(self.N_measurement):
            print(i)
            self.sweep()
            measurements[i] = observable(self.lattice)
        results = [0 for i in range(11)]
        for step in steps:
            result = 0
            for i in range(self.N_measurement-step):
                result += measurements[i] * measurements[i+step]
            result = (result/(self.N_measurement-step) -np.average(measurements)**2)/ (np.average([m**2 for m in measurements])-np.average(measurements)**2)
            results[step] = result 
        return results
        
def main():
    lat = Ising(N=20,kT=100000,N_thermal=100,N_measuremnt=2*1000,N_sweeps=1)
    print(lat.autocorrelation(Ising.get_magnetisation))
    
main()