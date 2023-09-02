import numpy as np
import matplotlib.pyplot as plt
from Stats import Stats

def analytical(a,N,mu,m):

    R = 1+ a**2*mu**2/(2*m)-a*mu/np.sqrt(m)*(1+a**2*mu**2/(4*m))**0.5

    return 1/(2*mu*(m+a**2*mu**2/4)**0.5)*(1+R**N)/(1-R**N)

class Lattice(object):

    def __init__(self, N, a, N_thermal, N_measurement, N_sweeps, Delta, N_tau, d_tau, mass, w) -> None:

        self.m = mass 

        self.w = w

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
    def action(lattice,a, m, w):

        S = 0

        for i in range(len(lattice)):

            S+= (1/2*a*( (lattice[(i+1)%len(lattice)] -lattice[i])/a)**2 +1/2 *a*m*w**2*lattice[i]**2)

        return S
    
    def HMC(self):

        ps = [np.random.normal() for i in range(self.N)]

        H = sum([p**2 for p in ps])/2 + Lattice.action(self.lattice,self.a,self.m,self.w)

        p_f, x_f = self.molecular_dynamics(ps.copy(),self.lattice.copy())

        H_new = sum([p**2 for p in p_f])/2 + Lattice.action(x_f,self.a,self.m,self.w)

        delta_H = H_new- H

        if delta_H <0 or np.exp(-delta_H) > np.random.uniform(0,1):
            self.lattice = x_f.copy()

    def molecular_dynamics(self,p_0, x_0):

        p = [p_0[i] - self.d_tau/2 * (self.m*1/self.a * (2*x_0[i]-x_0[(i-1)%self.N] - x_0[(i+1)%self.N]) + self.a*self.m*self.w**2 *x_0[i])  for i in range(self.N)]

        x = [x_0[i] + self.d_tau* p[i] for i in range(self.N)]

        for j in range(self.N_tau):

            p = [p[i] - self.d_tau * (self.m*1/self.a  * (2*x[i]-x[(i-1)%self.N] - x[(i+1)%self.N]) + self.a*self.m*self.w**2 *x[i])  for i in range(self.N)]

            x = [x[i] + self.d_tau* p[i] for i in range(self.N)]

        p = [p[i] - self.d_tau/2 * (self.m*1/self.a * (2*x[i]-x[(i-1)%self.N] - x[(i+1)%self.N]) + self.a*self.m*self.w**2 *x[i])  for i in range(self.N)]

        return p, x
    
    def thermalize(self):

        for i in range(self.N_thermal):
            print(i)

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

            if i% 1000 == 0:
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
        X = 20
        values = [0 for l in range(X)]

        
        for n in range(X):
            values[n] = 1/N* sum([lattice[((j+n)%N)]* lattice[j] for j in range(N)])

        

        return values

    @staticmethod
    def measure_position(lattice):
        return np.average(lattice)
    @staticmethod
    def measure_sq_position(lattice):
        return np.average([x**2 for x in lattice])
    
    @staticmethod
    def save_measurements(N, a, N_thermal, N_meausre, N_sweeps, measurements, observable_name, file_name):
        file = open(file_name,'a')
        
        file.write("########################"+ '\n')
        
        file.write("Description:" + '\n')

        file.write('The Data file that comes with this one stores the measremts of the observable ' + observable_name+ ". The following lines give a breif desription of the simulation used to draw them.")
        
        file.write('The lattice has '+ str(N)+ "sites and a is chosen to be " + str(a) + ' in this simulation'+ '\n')
        
        file.write(str(N_thermal) + ' sweeps have been at the beginning to thermalize the lattice'+ '\n')
        
        file.write(str(N_meausre) + ' Configurations of the lattice have been saved'+ '\n')
        
        file.write(str(N_sweeps)+ ' sweeps have been performed between each saved configuration to minimize the autocorrelation'+ '\n')
        
        file.write("########################"+ '\n')
        
        file.close()

        np.save(file_name,measurements)

def main():
    
   
    file_name = 'Test 3 B QSHO x^2'
    observalbe_name = 'x^2'
    a =  0.4
    N = 100
    
    N_measurements = 10**5

    N_sweeps = 1

    N_thermal =100

    Delta = 1

    w = 1
    
    m = 1


    lat = Lattice(N = N, a = a, N_thermal = N_thermal, N_measurement = N_measurements, N_sweeps = N_sweeps, Delta = Delta,N_tau=100,d_tau=0.01, mass=m, w=w)

    values = lat.generate_measurements(Lattice.measure_greens)
    cutoff = 15
    Stat = Stats(values)
    #print(Stat.estimate(cutoff=cutoff))
    X = 20

    vals = [[values[i][n] for i in range(N_measurements)] for n in range(20)]
    Stat = [Stats(b) for b in vals]
    estimates = [S.estimate(cutoff) for S in Stat]
    print(estimates)
    value_estimates = [estimates[i][0] for i in range(X)]
    error_estimates = [estimates[i][1] for i in range(X)]
    plt.errorbar([i for i in range(X)],value_estimates,yerr=error_estimates,fmt='.k', label = 'Simulation Data')
    plt.show()

    

    """ 
    result = [0 for i in range(t)]
    t = 10
    poss_a = np.linspace(0.1,1,t)
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
    
    plt.show()"""
    

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
