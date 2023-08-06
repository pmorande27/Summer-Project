import numpy as np
import matplotlib.pyplot as plt
class Stats(object):

    def __init__(self, lattices) -> None:

        self.lattices = lattices

    
    def measure(self, observable):

        return [observable(config) for config in self.lattices]
    
    def autocorrelation(self, observable,plot= False):

        steps = [i for i in range(20)]

        measurements = self.measure(observable=observable)

        average = np.average(measurements)

        average_of_sq = np.average([m**2 for m in measurements])

        sigma_sq = (average_of_sq-average**2)

        results = [0 for i in range(20)]

        for step in steps:

            result = 0

            for i in range(len(self.lattices)-step):

                result += (measurements[i]-average)*(measurements[i+step]-average)
            
            result = (result/(len(self.lattices)-step))/ sigma_sq
            
            results[step] = result 
        
        if plot:

            plt.plot(results,'o')

            plt.show()
        
        return results
        
    def integrated_autoccorelation(self, observable):

        return 2*(1/2 + sum(self.autocorrelation(observable=observable)[1:]))
    
    def estimate(self, observable):

        measurements = self.measure(observable=observable)

        average = np.average(measurements)

        average_of_sq = np.average([m**2 for m in measurements])

        sigma_sq = (average_of_sq-average**2)

        integrated_autoccorelation = self.integrated_autoccorelation(observable=observable)

        true_variance = sigma_sq*integrated_autoccorelation

        error = true_variance**0.5

        return average, error


