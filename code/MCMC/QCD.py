import numpy as np
def main():
    l = Lattice(1)

class Lattice(object):

    def __init__(self,N) -> None:
        self.N_x = N
        self.N_y = N
        self.N_z = N
        self.N_t = N
        U = np.zeros((self.N_t,self.N_x,self.N_y,self.N_z,4,3,3),dtype=complex) 
        for t in range(N):
            for x in range(N):
                for y in range(N):
                    for z in range(N):
                        for u in range(4):
                            for n in range(3):
                                U[t,x,y,z,u,n,n] = 1
                            print(U[t,x,y,z,u])
main()