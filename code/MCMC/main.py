from Stats import Stats
import numpy as np
import matplotlib.pyplot as plt

def analytical(a,N,mu,m):

    R = 1+ a**2*mu**2/(2*m)-a*mu/np.sqrt(m)*(1+a**2*mu**2/(4*m))**0.5

    return 1/(2*mu*(m+a**2*mu**2/4)**0.5)*(1+R**N)/(1-R**N)
def main():
    values = []
    err = []
    a = [0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05]
    
    for i in range(1,11):
        file_name = 'Results/Test '+str(i) + ' QSHO x^2.npy'
        measurements = np.load(file_name)
        Statistics = Stats(measurements)
        val1,err1 = Statistics.estimate(15)
        values += [val1]
        err += [err1]
       
        file_name = 'Results/Test '+str(i) + ' B QSHO x^2.npy'
        measurements = np.load(file_name)
        Statistics = Stats(measurements)
        val2,err2 = Statistics.estimate(15)
        values += [val2]
        err += [err2]
    plt.errorbar(a,values,yerr=err,fmt='.k', label = 'Simulation Data')
    
    pos_a = np.linspace(0.1,1.05,1000)
    real =[ analytical(a,100,1,1) for a in pos_a ]
    plt.plot(pos_a,real,label = "Theory")
    plt.xlabel('a')
    plt.ylabel('<x^2>')
    plt.legend()
    plt.show()
    
 

"""    
def generate_data():

    N = 8

    beta = 5.5

    N_matrix = 100

    epsilon = 0.24

    hits = 10

    N_measure = 100

    N_corr = 2

    N_thermal = 100

    lat = Lattice(N=N, beta = beta, Nmatrix= N_matrix, epsilon= epsilon, hits= hits, N_measruements=N_measure, N_corr= N_corr, N_themral=N_thermal)

    measurements = lat.generate_measurements(Lattice.average_plaquette)

    Lattice.save_measurements(N,beta,N_thermal,N_measure,N_corr,hits,epsilon,N_matrix, measurements, 'Average Plaquette', 'Measurement 1')
"""


main()


