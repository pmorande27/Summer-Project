import numpy as np
import random
import cmath
import math

def main():
    l = Lattice(8,6,5,0.24,10)

    for i in range(10):
        l.sweep()
        print(i)

class Lattice(object):

    def __init__(self,N,beta,Nmatrix,epsilon, hits) -> None:
        self.N = N
        self.N_x = N
        self.N_y = N
        self.N_z = N
        self.N_t = N
        self.Nmatrix = Nmatrix
        self.beta = beta
        self.epsilon = epsilon
        self.hits = hits
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
    def dagger(self,D):
        N=len(D)
        H=np.zeros((N,N),dtype=complex)
        R=np.matrix(D)
        H=R.getH()
        return H.copy()
    

    def generte_random_matrix(self):
        identity=np.eye(3)
        H=np.zeros((3,3),dtype=complex)
        w=cmath.sqrt(-1)
        for s in range(self.Nmatrix):
            for j in range(0,3):
                for i in range(0,3):
                    H[j,i]=complex(random.uniform(-1,1),random.uniform(-1,1))
            H=(H.copy()+self.dagger(H.copy()))/2. #generation o the hermitian matrice
            for n in range(30): #Taylor series for a SU(3) mstrice
                self.M[s]=self.M[s]+(w*self.epsilon)**n/math.factorial(n)*np.linalg.matrix_power(H,n)
            self.M[s]=self.M[s]/np.linalg.det(self.M[s])**(1/3) #Unitarization

            self.M[s+self.Nmatrix]=self.dagger(self.M[s]) #Saving the inverse

    def calculate_gamma(self,position_link):
        result = np.zeros((3,3),dtype=complex)
        t,x,y,z,u = position_link[0],position_link[1],position_link[2],position_link[3],position_link[4]
        for v in range(4):
            if v != u:
                vs = [0 for i in range(4)]
                vs[v] = 1
                us = [0 for i in range(4)]
                us[u] = 1
                part_one = np.dot(np.dot(self.U[(t+us[0])%self.N,(x+us[1])%self.N,(y+us[2])%self.N,(z+us[3])%self.N,v],self.dagger(self.U[(t+vs[0])%self.N,(x+vs[1])%self.N,(y+vs[2])%self.N,(z+vs[3])%self.N,u])),self.dagger(self.U[t,x,y,z,v]))
                part_two = np.dot(np.dot(self.dagger(self.U[(t+us[0]-vs[0])%self.N,(x+us[1]-vs[1])%self.N,(y+us[2]-vs[2])%self.N,(z+us[3]-vs[3])%self.N,v]),self.dagger(self.U[(t-vs[0])%self.N,(x-vs[1])%self.N,(y-vs[2])%self.N,(z-vs[3])%self.N,u])),self.U[(t-vs[0])%self.N,(x-vs[1])%self.N,(y-vs[2])%self.N,(z-vs[3])%self.N,v])
                result += part_one + part_two
        return result
    def sweep(self):
        for t in range(self.N_t):
            for x in range(self.N_x):
                for y in range(self.N_y):
                    for z in range(self.N_z):
                        for u in range(4):
                            position_link = np.array([t,x,y,z,u])
                            Gamma = self.calculate_gamma(position_link)
                            for hit in range(self.hits):
                                s = random.randint(0,2*self.Nmatrix-1)
                                M =self.M[s]
                                Delta_U = np.dot(M,self.U[t,z,y,z,u])
                                Delta_S =  -self.beta/3 * np.real(np.trace(np.dot((Delta_U-self.U[t,z,y,z,u]),Gamma)))
                                prob = min(1,np.exp(-Delta_S))
                                if random.random() <= prob:
                                    self.U[t,x,y,z,u] = Delta_U

    def get_plaquette(self,t,x,y,z,u,v):
        pass


    def action(self):
        result =0
        for t in range(self.N_t):
            for x in range(self.N_x):
                for y in range(self.N_y):
                    for z in range(self.N_z):
                        for v in range(4):
                            for u in range(v+1,4):
                                result = self.beta*(1-self.get_plaquette(t,x,y,z,u,v))



main()