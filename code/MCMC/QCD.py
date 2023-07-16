import numpy as np
import random
def main():
    l = Lattice(1)

class Lattice(object):

    def __init__(self,N,beta) -> None:
        self.N_x = N
        self.N_y = N
        self.N_z = N
        self.N_t = N
        self.beta = beta
        self.U = np.zeros((self.N_t,self.N_x,self.N_y,self.N_z,4,3,3),dtype=complex) 
        for t in range(N):
            for x in range(N):
                for y in range(N):
                    for z in range(N):
                        for u in range(4):
                            for n in range(3):
                                self.U[t,x,y,z,u,n,n] = 1
                            print(self.U[t,x,y,z,u])
    
    def generte_random_matrix(self):
        pass
    def calculate_gamma(self,position_link):
        pass

    def sweep(self):
        for t in range(self.N_t):
            for x in range(self.N_x):
                for y in range(self.N_y):
                    for z in range(self.N_z):
                        for u in range(4):
                            position_link = np.array([t,x,y,z,u])
                            M = self.generte_random_matrix()
                            Delta_U = np.dot(M,self.U[t,z,y,z,u])
                            Gamma = self.calculate_gamma(position_link)
                            Delta_S =  -self.beta/3 * np.real(np.trace(np.dot((Delta_U-self.U[t,z,y,z,u]),Gamma)))
                            prob = min(1,np.exp(-Delta_S))
                            if random.random() <= prob:
                                self.U[t,x,y,z,u] = Delta_U
                                



main()