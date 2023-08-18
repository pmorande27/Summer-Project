import numpy as np
import matplotlib.pyplot as plt

def analytical(a,N,mu,m):

    R = 1+ a**2*mu**2/(2*m)-a*mu/np.sqrt(m)*(1+a**2*mu**2/(4*m))**0.5

    return 1/(2*mu*(m+a**2*mu**2/4)**0.5)*(1+R**N)/(1-R**N)

class Lattice(object):

    def __init__(self, N, a, N_thermal, N_measurement, N_sweeps, Delta, N_tau, d_tau) -> None:

        self.N_thermal = N_thermal
        
        self.N_sweeps = N_sweeps
        
        self.N_measurement = N_measurement
        
        self.a = a
        
        self.N = N
        
        self.Delta = Delta
        
        self.lattice = np.zeros(N)
        
        self.total_time = N*a
        
        self.sweeps = 0
        
        self.N_tau = N_tau
        
        self.d_tau = d_tau
        
        self.thermalize()
    
    @staticmethod
    def action(lattice,a):

        m = 1

        mu = 1

        S = 0

        for i in range(len(lattice)):

            S+= a*(1/2*( (lattice[(i+1)%len(lattice)] -lattice[i])/a)**2 +1/2 *lattice[i]**2)

        return S
    
    def HMC(self):

        ps = [np.random.normal() for i in range(self.N)]

        H = sum([p**2 for p in ps])/2 + Lattice.action(self.lattice,self.a)

        p_f, x_f = self.molecular_dynamics(ps.copy(),self.lattice.copy())

        H_new = sum([p**2 for p in p_f])/2 + Lattice.action(x_f,self.a)

        delta_H = H_new- H

        if delta_H <0 or np.exp(-delta_H) > np.random.uniform(0,1):
            self.lattice = x_f.copy()

    def molecular_dynamics(self,p_0, x_0):

        p = [p_0[i] - self.d_tau/2 * (1/self.a * (2*x_0[i]-x_0[(i-1)%self.N] - x_0[(i+1)%self.N]) + self.a *x_0[i])  for i in range(self.N)]

        x = [x_0[i] + self.d_tau* p[i] for i in range(self.N)]

        for j in range(self.N_tau):

            p = [p[i] - self.d_tau * (1/self.a * (2*x[i]-x[(i-1)%self.N] - x[(i+1)%self.N]) + self.a *x[i])  for i in range(self.N)]

            x = [x[i] + self.d_tau* p[i] for i in range(self.N)]

        p = [p[i] - self.d_tau/2 * (1/self.a * (2*x[i]-x[(i-1)%self.N] - x[(i+1)%self.N]) + self.a *x[i])  for i in range(self.N)]

        return p, x
    
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
    
    def generate_measurements(self, observable):

        results = [0 for i in range(self.N_measurement)]

        for i in range(self.N_measurement):

            for j in range(self.N_sweeps):

                self.HMC()

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
    
    t = 10
    
    poss_a = np.linspace(0.1,1,t)
    
    result = [0 for i in range(t)]
    
    for i in range(len(poss_a)):
    
        N = 100
    
        a = poss_a[i]
    
        N_measurements = 10**4
    
        N_sweeps = 1
    
        N_thermal =100
    
        Delta = 1
    
        lat = Lattice(N = N, a = a, N_thermal = N_thermal, N_measurement = N_measurements, N_sweeps = N_sweeps, Delta = Delta,N_tau=100,d_tau=0.01)
    
        values = lat.generate_measurements(Lattice.measure_sq_position)
    
        result[i] = np.average(values)
    
    real = [analytical(a,100,1,1) for a in np.linspace(0.1,1,100)]
    
    plt.plot(poss_a,result,'o')
    
    plt.plot(np.linspace(0.1,1,100),real)
    
    plt.show()
    

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

