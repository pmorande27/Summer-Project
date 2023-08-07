import numpy as np
import cmath
import math
from scipy.stats import sem


def main():
    Nmeasure = 10
    Ncorr = 25
    N = 8
    beta = 5.5
    N_matrix = 100
    epsilon = 0.24
    hits = 10
    N_thermal = 100
    
    print("Number of Measurments: "+ str(Nmeasure))
    print("Number of sweeps before Correlation: " + str(Ncorr))

    l = Lattice(N, beta, N_matrix, epsilon, hits, Nmeasure, Ncorr, N_thermal)
    avg = l.generate_configurations()
    Lattice.save_results(N,beta,N_thermal,Nmeasure,Ncorr,hits,epsilon,N_matrix,avg, 'Test 2')
    

class Lattice(object):

    def __init__(self,N,beta,Nmatrix,epsilon, hits, N_measruements, N_corr,N_themral) -> None:

        self.N,self.N_x,self.N_y,self.N_z,self.N_t = N, N, N, N, N

        

        self.Nmatrix, self.beta, self.epsilon, self.hits, self.N_thermal = Nmatrix, beta, epsilon, hits, N_themral

        self.N_measurements, self.N_corr = N_measruements, N_corr
    
        self.M=np.zeros((Nmatrix*2,3,3),dtype=complex)
        
        self.generte_random_matrix()

        self.U = np.zeros((self.N_t,self.N_x,self.N_y,self.N_z,4,3,3),dtype=complex) 
        

        for t in range(N):
        
            for x in range(N):
        
                for y in range(N):
        
                    for z in range(N):
        
                        for u in range(4):
        
                            for n in range(3):
        
                                self.U[t,x,y,z,u,n,n] = 1
        
        self.thermalize()

    

    def generte_random_matrix(self):

        identity=np.eye(3)

        H=np.zeros((3,3),dtype=complex)

        w=cmath.sqrt(-1)

        for s in range(self.Nmatrix):

            for j in range(0,3):

                for i in range(0,3):

                    H[j,i]=complex(np.random.uniform(-1, 1), np.random.uniform(-1, 1))

            H=(H.copy()+Lattice.dagger(H.copy()))/2. #generation o the hermitian matrice

            for n in range(30): #Taylor series for a SU(3) mstrice

                self.M[s] = self.M[s] + (w*self.epsilon)**n/math.factorial(n)*np.linalg.matrix_power(H, n)

            self.M[s] = self.M[s]/np.linalg.det(self.M[s])**(1/3) #Unitarization

            self.M[s + self.Nmatrix] = Lattice.dagger(self.M[s]) #Saving the inverse

    def calculate_gamma(self, position_link):

        t, x, y, z, ro, N = position_link[0], position_link[1], position_link[2], position_link[3], position_link[4], self.N

        gamma = 0. #inizializing gamma
        
        inc_ro = np.zeros((4), 'int')
        
        inc_ni = np.zeros((4), 'int')
        
        position = np.array([t, x, y, z])
        
        inc_ro[ro] = 1 #increment in the ro direction
        #next site on mi
        position_one = (position + inc_ro) % N

        for ni in range(0, 4):

            if ni != ro :

                inc_ni[ni] = 1 #increment in the ni direction
                #next site on ni
                position_two = (position + inc_ni) % N
                #Previous site on ni
                position_four = (position - inc_ni) % N
                #next site on mi and previous on ni
                position_three = (position + inc_ro - inc_ni) % N

                inc_ni[ni] = 0

                gamma+= np.dot(np.dot(self.U[tuple(position_one)+(ni,)],Lattice.dagger(self.U[tuple(position_two)+(ro,)])),Lattice.dagger(self.U[t,x,y,z,ni])) \
                        \
                        +np.dot(np.dot(Lattice.dagger(self.U[tuple(position_three)+(ni,)]),Lattice.dagger(self.U[tuple(position_four)+(ro,)])),self.U[tuple(position_four)+(ni,)])
       
        return gamma.copy()
    
    def sweep(self):

        for t in range(self.N_t):

            for x in range(self.N_x):

                for y in range(self.N_y):

                    for z in range(self.N_z):

                        for u in range(4):

                            position_link = np.array([t,x,y,z,u])

                            Gamma = self.calculate_gamma(position_link)

                            for hit in range(self.hits):

                                s = np.random.randint(0,2*self.Nmatrix)

                                M =self.M[s].copy()

                                Delta_U = np.dot(self.M[s].copy(),self.U[t,x,y,z,u].copy())

                                Delta_S = -self.beta/(3)*np.real(np.trace(np.dot((np.dot(self.M[s].copy(),self.U[t,x,y,z,u].copy())-self.U[t,x,y,z,u].copy()),Gamma.copy())))

                                if Delta_S <0 or np.exp(-Delta_S) > np.random.uniform(0,1):

                                    self.U[t,x,y,z,u] = Delta_U.copy()

    def thermalize(self):

        for i in range(self.N_thermal):

            self.sweep()

            print(i)

            print(Lattice.average_plaquette(self.U))

    def generate_configurations(self):
        results = [0 for i in range(self.N_measurements)]
        for i in range(self.N_measurements):
            
            for j in range(self.N_corr):
                self.sweep()
                print(j)
            results[i] = self.U.copy()
        return(results)
    @staticmethod
    def get_plaquette(t,x,y,z, U):

        N, WL =len(U), 0

        

        incmi, incni =np.zeros((4),'int'), np.zeros((4),'int')

        position = np.array([t, x, y, z])

        for mi in range(0, 4):

            incmi[mi] = 1 #increment in the mi direction
        #next site on mi
            pos_one = (position + incmi) %N
            for ni in range(0,mi):
                incni[ni]=1 #increment in the ni direction
                #next site on ni
                pos_two = (position + incni) %N 

                incni[ni]=0
            
                WL+=np.trace(np.dot(U[tuple(position)+(mi,)],np.dot(np.dot(U[tuple(pos_one)+(ni,)],Lattice.dagger(U[tuple(pos_two)+(mi,)])),Lattice.dagger(U[tuple(position)+(ni,)]))))
            
            incmi[mi]=0

        return np.real(WL)/(3.*6.)


    @staticmethod
    def average_plaquette(lattice):
        N = len(lattice)
        result =0

        for t in range(N):

            for x in range(N):

                for y in range(N):

                    for z in range(N):
                            
                            result += Lattice.get_plaquette(t, x, y, z, lattice)

        return (result)/(N**4)

    @staticmethod
    def dagger(D):

        N = len(D)

        H = np.zeros((N, N), dtype=complex)

        R=np.matrix(D)

        H=R.getH()

        return H.copy()
    
    @staticmethod
    def save_results(N, beta, N_thermal, N_meausre, N_correlation, hits, epsilon, N_matrix, lattices, file_name):
        file = open(file_name,'a')
        file.write("########################"+ '\n')
        file.write("Description:" + '\n')
        file.write('The lattice has '+ str(N)+'**4 sites and beta is chosen to be ' + str(beta) + ' in this simulation'+ '\n')
        file.write(str(N_thermal) + ' sweeps have been at the beginning to thermalize the lattice'+ '\n')
        file.write(str(N_meausre) + ' Configurations of the lattice have been saved'+ '\n')
        file.write(str(N_correlation)+ ' sweeps have been performed between each saved configuration to minimize the autocorrelation'+ '\n')
        file.write('Each link is updated '+ str(hits) +' times each time that it is selected'+ '\n')
        file.write('The parameter that controls the deviation form the identity matrix of the random SU(3) elements has been ' + str(epsilon) + ' in this simulation and ' + str(N_matrix) + ' random SU(3) matrices have been generated'+ '\n')
        file.write("########################"+ '\n')
        file.close()
        np.save(file_name,lattices)
     
    
    
    """def obtain_spatial_loops(self, t):
        WL =0

        for x in range(self.N_x):

            for y in range(self.N_y):

                for z in range(self.N_z):
                     
                    N =self.N 

                    incmi, incni =np.zeros((4),'int'), np.zeros((4),'int')
                     
                    position = np.array([t, x, y, z])
                    for mi in range(1, 4):

                        incmi[mi] = 1 
                        pos_one = (position + incmi) %N
                        for ni in range(1,mi):
                            incni[ni]=1 
                   
                            pos_two = (position + incni) %N 

                            incni[ni]=0
                
                            WL+=np.trace(np.dot(self.U[tuple(position)+(mi,)],np.dot(np.dot(self.U[tuple(pos_one)+(ni,)],Lattice.dagger(self.U[tuple(pos_two)+(mi,)])),Lattice.dagger(self.U[tuple(position)+(ni,)]))))
                
                        incmi[mi]=0



        return (np.real(WL)/(3.*3))/self.N**3
    

    def measure_at_diff_times(self):

        results = [0 for i in range(self.N_t)]

        for t in range(self.N_t):
            results[t] = self.obtain_spatial_loops(t)
        
        return results
    

    def measure_two(self):
        results = [0 for i in range(self.N_measurements)]
        for i in range(self.N_measurements):
            
            for j in range(self.N_corr):
                self.sweep()
                print(j)
            results[i] = self.measure_at_diff_times()
        return(results)"""
#main()