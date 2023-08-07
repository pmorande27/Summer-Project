from Stats import Stats
import numpy as np
from QCD import Lattice
def main():
    lattices = np.load('Test 2.npy')
    Statistics = Stats(lattices)
    print(Statistics.estimate(Lattice.average_plaquette,3))
main()
