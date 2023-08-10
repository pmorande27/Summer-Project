import numpy as np
import matplotlib.pyplot as plt
class Stats(object):

    def __init__(self, measurements) -> None:

        self.measurements = measurements

    
    def autocorrelation(self, cutoff, plot= False):

        steps = [i for i in range(cutoff)]

        measurements = self.measurements

        average = np.average(measurements)

        sigma_sq = np.var(measurements,ddof=1)

        #sigma_sq = (average_of_sq-average**2)

        results = [0 for i in range(cutoff)]

        for step in steps:

            result = 0

            for i in range(len(self.measurements)-step):

                result += (measurements[i]-average)*(measurements[i+step]-average)
            
            result = (result/(len(self.measurements)-step))/ sigma_sq
            
            results[step] = result 
        
        if plot:

            plt.plot(results,'o')

            plt.show()
        
        return results
        
    def integrated_autoccorelation(self, cutoff):
        a = self.autocorrelation(cutoff=cutoff)[1:]

        return 2*(1/2 + sum(self.autocorrelation(cutoff=cutoff)[1:]))
    
    def estimate(self, cutoff):

        measurements = self.measurements
        average = np.average(measurements)

        average_of_sq = np.average([m**2 for m in measurements])

        sigma_sq = np.var(measurements,ddof=1)

        integrated_autoccorelation = self.integrated_autoccorelation(cutoff=cutoff)

        true_variance = sigma_sq*integrated_autoccorelation

        error = np.sqrt(true_variance/len(measurements))

        print(integrated_autoccorelation)

        return average, error
    @staticmethod
    def get_measurements(file_name, observable_name, observable):
        lattices = np.load(file_name)
        name = file_name.split('.')[0] + " " + observable_name
        results = [observable(config) for config in lattices]
        np.savetxt(name,results)
