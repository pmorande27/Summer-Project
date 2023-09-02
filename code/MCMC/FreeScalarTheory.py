import numpy as np

class Lattice(object):

    def __init__(self, N, k, l, a, N_measurement, N_sweeps, N_thermal, N_tau, d_tau) -> None:
        
        self.N = N
        self.lattice = np.zeros((N,N,N))
        self.N_measurement = N_measurement
        self.N_sweeps = N_sweeps
        self.N_thermal = N_thermal
        self.k = k
        self.l = l
        self.a = a
        self.N_tau = N_tau
        self.d_tau = d_tau
        self.themralize()
    @staticmethod
    def action(lattice, k, l):
        action = 0
        N = len(lattice)
        for t in range(N):
            for x in range(N):
                for y in range(N):
                        action += -2*k *lattice[t,x,y]* ( lattice[t,x,(y+1)%N]+  lattice[t,(x+1)%N,y]+ lattice[(t+1)%N,x,y])
                        action += lattice[t,x,y]**2 + l*(lattice[t,x,y]**2-1)**2
        return action
    @staticmethod
    def calculate_hamiltonian(p,phi,k,l):
        N = len(p)
        H = 0
        for t in range(N):
            for x in range(N):
                for y in range(N):
                    H += 1/2 * p[t,x,y]**2
        H += Lattice.action(phi, k, l)
        return H
    
    @staticmethod
    def move_p(p,phi,eps,k,l):
        N = len(p)
        for t in range(N):
            for x in range(N):
                for y in range(N):
                        p[t,x,y] = p[t,x,y]- Lattice.force(phi,t,x,y,k,l)*eps
    
    @staticmethod
    def move_phi(p, phi, eps):
        N = len(p)
        for t in range(N):
            for x in range(N):
                for y in range(N):
                        phi[t,x,y] = phi[t,x,y]+ p[t,x,y]*eps

    @staticmethod
    def force(phi,t,x,y, k , l):
        N = len(phi)
        force = 0
        force +=  -2*k *(  phi[t,x,(y+1)%N]+  phi[t,(x+1)%N,y]+ phi[(t+1)%N,x,y] + phi[t,x,(y-1)%N]+  phi[t,(x-1)%N,y]+ phi[(t-1)%N,x,y])
        force += 2*phi[t,x,y]
        force += 4*l*(phi[t,x,y]**2-1)*phi[t,x,y]
        return force
    
    def measure_magnetization(lattice):
        return np.sum(lattice)
    def HMC(self):
        p = np.random.normal(size=(self.N,self.N,self.N))

        phi = self.lattice.copy()

        H = Lattice.calculate_hamiltonian(p,phi.copy(),self.k,self.l)
        

        p_f, phi_f = self.molecular_dynamics(p.copy(),phi)

        H_new = Lattice.calculate_hamiltonian(p_f,phi_f,self.k,self.l)
        delta_H = H_new- H

        if delta_H <0 or np.exp(-delta_H) > np.random.uniform(0,1):
            self.lattice = phi_f.copy()




    
    def molecular_dynamics(self, p, phi):
        Lattice.move_p(p,phi,self.d_tau/2,self.k,self.l)
        for i in range(1,self.N_tau):
            Lattice.move_phi(p,phi,self.d_tau)
            Lattice.move_p(p,phi,self.d_tau/2,self.k,self.l)

        Lattice.move_phi(p,phi,self.d_tau)
        Lattice.move_p(p,phi,self.d_tau/2,self.k,self.l)



        return p, phi



    def themralize(self):
        for i in range(self.N_thermal):
            print(i)
            print(Lattice.measure_magnetization(self.lattice))
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
    
def main():
    lat = Lattice(6,0.185825, 1.1689,0.1, 10, 1, 100, 1000, 0.01)
    b = lat.generate_measurements(Lattice.measure_magnetization)
    print(np.average(b))

main()