import numpy as np
import matplotlib.pyplot as plt
class Gaussian(object):
    def __init__(self,epsilon, N_thermal, N_measurment, N_sweeps, N_tau, d_tau) -> None:

        self.lattice = 0

        self.N_thermal = N_thermal
        
        self.N_measurement = N_measurment
        
        self.N_sweeps = N_sweeps
        
        self.epsilon = epsilon

        self.N_tau = N_tau
        
        self.d_tau = d_tau
        
        self.thermalize()
    def molecular_dynamics(self, p_0, x_0):
        x = x_0
        p = p_0
        x = x + p *self.d_tau/2
        for i in range(1,self.N_tau):
            p = p - x*self.d_tau
            x = x+ p*self.d_tau
        p = p - x*self.d_tau
        x = x + p*self.d_tau/2
        return p,x
        
    def HMC(self):
        p = np.random.normal()
        H = 0.5*p**2 + Gaussian.s(self.lattice)

        p_f, x_f = self.molecular_dynamics(p, self.lattice)
        H_new = 0.5*p_f**2 + Gaussian.s(x_f)
        delta_H = H_new- H
        if delta_H <0 or np.exp(-delta_H) > np.random.uniform(0,1):
    
            self.lattice = x_f


    @staticmethod
    def s(x):
    
        return x**2/2
    

    def thermalize(self):
    
        for i in range(self.N_thermal):
    
            self.HMC()
    
    def Get_configurations(self):
    
        results = [0 for i in range(self.N_measurement)]
    
        for i in range(self.N_measurement):
    
            for j in range(self.N_sweeps):
    
                self.HMC()
        
            results[i] = self.lattice
    
        return(results)


def main():

    Gauss = Gaussian(epsilon=1,N_thermal=100,N_measurment=10**7,N_sweeps=1,N_tau=40,d_tau=1)
    a = Gauss.Get_configurations()
    print(np.average(a))
    plt.hist(a,10000)
    plt.show()
    
    


main()
