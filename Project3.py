# Xavier Linn
# MSE 215, Project 3
# A Monte Carlo Problem
# Metropolis algorithm to compute the properties of a 2D Ising models.

# Import statements
import numpy as np
import matplotlib.pyplot as plt

# Class for creating a 2D array
class TwoDArray:

    # Instantiate a 2D array
    def __init__(self, size):
        self.size = size # 2D array is a square with length size. 0, 1, ... size - 1
        self.alloy = np.arange(size * size).reshape(size, size) # Generate a 2D Array
        self.populateTwoDArray() # Populate the array randomly with -1's and 1's

    # Populate alloy randomly with 1's an -1's
    def populateTwoDArray(self):
        # Iterate through all the 2D array and put a -1 or 1 based on random number
        # Generate a random number between 1 and 10, if < 6 populate with -1
        for i in range(0, self.size):
            for j in range(0, self.size):
                spin = 1
                randomValue = np.random.randint(1,11)
                if (randomValue <= 5):
                    spin = -1
                self.alloy[i][j] = spin

    # Prints the 2D array
    def printArray(self):
        print("Printing the 2D array")
        print(self.alloy)


    # Retrieve the state at index i
    def stateAtIndexI(self, i, j):
        if (i == -1):
            return self.alloy[self.size - 1][j]
        elif (i == self.size):
            return self.alloy[0][j]
        else:
            return self.alloy[i][j]

    # Retrieve the state at index j
    def stateAtIndexJ(self, i, j):
        if (j == -1):
            return self.alloy[i][self.size - 1]
        elif (j == self.size):
            return self.alloy[i][0]
        else:
            return self.alloy[i][j]

    # Change state at index i, j
    def changeState(self, site):
        i = site[0]
        j = site[1]
        currentState = self.alloy[i][j]
        if currentState == 1:
            self.alloy[i][j] = -1
        else:
            self.alloy[i][j] = 1

# Class for running the MonteCaro Simulation, one class instance
# per temp and system size
class MonteCarlo:

    # Instantiate a  Monte Carlo Simulation
    def __init__(self, s, temp):
        self.twoDArray = TwoDArray(s)
        self.totalEnergy = 0
        self.size = s
        self.time = 0
        self.temperature = temp
        self.boltzmanConstant = 1
        self.beta = 1.0 / (self.temperature * self.boltzmanConstant)
        self.avgMagArray = np.array
        self.magneticMomentArray = np.array
        self.energyArray = np.array
        self.steps = 0
        self.ntrials = 0
        self.correlationTime = 0
        self.correlationFxn = np.array

    # Plot total energy vs time
    def totalEnergyVsTime(self):
        plt.plot(self.energyArray)
        plt.title("ENERGY VS TIME")
        plt.show()

    # Plot avg mag vs time
    def magVsTime(self):
        plt.plot(self.avgMagArray)
        plt.title("MAG VS TIME")
        plt.show()

    # Calculate the initial total energy, i.e. before any spins have been flipped
    def calculateInitialTotalEnergy(self):
        for i in range(0, self.size):
            for j in range(0, self.size):
                self.totalEnergy += self.twoDArray.alloy[i][j]*(self.twoDArray.stateAtIndexI(i + 1, j) + self.twoDArray.stateAtIndexJ(i, j + 1))
        self.totalEnergy *= (-1)

    # Returns a random site as a tuple (i, j)
    def generateRandomSite(self):
        i = np.random.randint(0, self.size)
        j = np.random.randint(0, self.size)
        return (i, j)

    # Calculate change in energy:
    # calculates the energy between two iterations
    def energyContributionAtSite(self, site):
        i = site[0]
        j = site[1]
        energyContribution = (-1) * self.twoDArray.alloy[i][j] * (self.twoDArray.stateAtIndexI(i + 1, j) + self.twoDArray.stateAtIndexJ(i, j + 1) + self.twoDArray.stateAtIndexI(i - 1, j) + self.twoDArray.stateAtIndexJ(i, j - 1))
        return energyContribution

    # Determines the energy after a Monte Carlo time step (a spin flip attemp)
    def calculateTotalEnergyAfterSpinFlip(self, site):
        initialEnergy = self.totalEnergy
        energyContributionAtSiteBeforeFlip = self.energyContributionAtSite(site)
        self.twoDArray.changeState(site)
        energyContributionAtSiteAfterFlip = self.energyContributionAtSite(site)
        finalEnergy = initialEnergy - energyContributionAtSiteBeforeFlip + energyContributionAtSiteAfterFlip
        return finalEnergy

    # Accept the spin flip or not
    def acceptFlipOrNot(self, newEnergy):
        energyChange = newEnergy - self.totalEnergy
        if (energyChange <= 0):
            return True
        else:
            boltzman = np.exp(-1 * self.beta * energyChange)
            eta = np.random.random_sample()
            if (eta < boltzman): # If eta < boltzman factor, ACCEPT
                return True
            else: # If eta >= boltzman factor, REJECT
                return False

    # Calculates the average magnetization per site
    def calculateAvgMagnetization(self):
        avg = abs((np.sum(self.twoDArray.alloy, None, float) + 0.0) / np.power(self.size, 2))
        return avg

    # Calculates the magnetic moment
    def calculateMagneticMoment(self):
        return np.sum(self.twoDArray.alloy, None, float)

    # Metropolis Alogrithm
    def metropolisAlgorithm(self, steps):
        print("Running MD Simulation")
        print("size: ", self.size)
        print("steps: ", steps)

        # Setup data structures
        numSpin = self.size * self.size # By convention one "time" step is L^2 spin flip attempts
        initialNumDataPoints = int(steps / numSpin)
        self.avgMagArray = np.arange(initialNumDataPoints, dtype=np.float)
        self.magneticMomentArray = np.arange(initialNumDataPoints)
        self.energyArray = np.arange(initialNumDataPoints)
        #self.avgMagArray = np.arange(250000, dtype=np.float)
        #self.magneticMomentArray = np.arange(250000)
        #self.energyArray = np.arange(250000)

        # Take the system to equilibrium and record data
        self.calculateInitialTotalEnergy() # Calculate the total energy before any spin flip attempts
        for i in range(0, steps):
            site = self.generateRandomSite() # Select a random site
            totalEnergyAfterSpinFlip = self.calculateTotalEnergyAfterSpinFlip(site) # Calculate the energy after a random spin is flipped
            # Accept or reject the spin flip
            flip = self.acceptFlipOrNot(totalEnergyAfterSpinFlip)
            if (flip):
                self.totalEnergy = totalEnergyAfterSpinFlip
            else:
                self.twoDArray.changeState(site)
            # Store energy and magnetization while system reaches equilibrium
            if (i % numSpin == 0):
            #if (i > 250000):
                if (i / numSpin == initialNumDataPoints):
                    break
                #print(self.calculateAvgMagnetization())
                self.avgMagArray[self.time] = self.calculateAvgMagnetization()
                self.magneticMomentArray[self.time] = self.calculateMagneticMoment()
                #print(self.totalEnergy)
                self.energyArray[self.time] = self.totalEnergy
                self.time += 1

        # Plot before removing non-equilibrium data
        #self.totalEnergyVsTime()
        #self.magVsTime()


        # Determine the correlation time.
        self.correlationFxn = self.estimated_autocorrelation(self.avgMagArray)
        self.correlationTime = self.integrateCorrelation(self.correlationFxn)
        print("Size of avgMagArray w/ non-eq. data: ", self.avgMagArray.size)
        print("Size of energyArray w/ non-eq. data: ", self.energyArray.size)
        # Remove non-equilibrium data using correlation time. "Collect Data" now that the system is at equilibrium
        #scaledCorrelationTime = int(25000 * self.correlationTime)
        scaledCorrelationTime = 2500
        #print("Correlation time; ", self.correlationTime)
        #print("Correlation time scaled; ", scaledCorrelationTime)
        self.avgMagArray = self.avgMagArray[scaledCorrelationTime:]
        self.energyArray = self.energyArray[scaledCorrelationTime:]
        self.magneticMomentArray = self.magneticMomentArray[scaledCorrelationTime:]
        print("Size of avgMagArray w/ only eq. data: ", self.avgMagArray.size)
        print("Size of energyArray w/ only eq. data: ", self.energyArray.size)

        # Determine the uncertainty
        #uncertInMag = self.bootStrappingMethod()
        #print("The uncertainty in avg mag is: ", uncertInMag)

        # Plot after removing non-equilibrium data
        #self.totalEnergyVsTime()
        #self.magVsTime()

    # Calculate the correlation data, outputs an array
    def estimated_autocorrelation(self, x):
        n = len(x)
        variance = x.var()
        x = x-x.mean()
        r = np.correlate(x, x, mode = 'full')[-n:]
        result = r/(variance*(np.arange(n, 0, -1)))
        return result

    # Approximate the index where the correlation function corsses the x-axis
    def findZeroIndex(self, x):
        i = 0
        j = 1
        for v in x: # Finds the two points where the value changes from (+) to (-)
            if x[i] > 0 and x[j] < 0:
                return j
            i += 1
            j += 1

    # Using the trapezoid method integrate the correlation function
    def integrateCorrelation(self, correlationFxn):
        index = self.findZeroIndex(correlationFxn)
        correlationFxn = correlationFxn[:index + 1]
        correlationTime = np.trapz(correlationFxn)
        return correlationTime

    # Calculates the number of independent trials using the correlation time
    def calculateNumberOfIndpTrials(self):
        #print(self.avgMagArray.size)
        self.ntrials = int(self.avgMagArray.size / (2 * self.correlationTime))

    # Bootstrapping method to determine the uncertainty in a given data set
    def bootStrappingMethod(self, data):
        self.calculateNumberOfIndpTrials()
        p = np.arange(self.ntrials, dtype='f')
        avgSampleArray = np.arange(1000, dtype='f')
        for i in range(0, 1000):
            for j in range(0, int(self.ntrials)):
                index = np.random.randint(1, data.size)
                p[j] = data[index]
            avgSample = p.mean()
            avgSampleArray[i] = avgSample
        avp = np.sqrt(np.var(avgSampleArray))
        return avp

    # Calculates the the variance
    def calculateHeatCapcaityPerSite(self):
        # Cv = (1 / (Kb * T^2 * L^2)) * (<E^2> - <E>^2)
        var = np.var(self.energyArray)
        Cv = (1.0 / (np.power(self.size, 2) * np.power(self.temperature, 2))) * var
        return Cv

    # Calculates the exact solution for the magnetic susceptibility
    def calculateExactMagneticSusceptibilityPerSite(self):
        avg_x = np.power((1 - np.power(np.sinh(2 * (1.0 / self.temperature) * 1), -4)), 1/8)
        return avg_x

    # Calculate the critical temperature
    def calcuateExactCriticalTemperature(self):
        return (2.0 / np.log(1.0 + np.sqrt(2)))

    # Calculates the exact heat capacity per spin
    def calculateExactHeatCapacity(self):
        Tc = self.calcuateExactCriticalTemperature()
        print("Tc: ", Tc)
        val = (2.0 / np.pi) * np.power((2.0 / Tc), 2) * (-1 * np.log(1 - (self.temperature + 0.0) / Tc) + np.log(Tc / 2.0) - (1 + (np.pi + 0.0)/ 4.0))
        return val

    # Caclulate the magnetic susceptibility per site
    def calculateMagneticSusceptibilityPerSite(self):
        # x = beta * (1 / L^2) (<m^2> - <m>^2), beta = 1 / T
        var = np.var(self.magneticMomentArray)
        print(var)
        #var = np.var(self.avgMagArray)
        x = self.beta * (1.0 / np.power(self.size, 2)) * var
        return x

def simulation():
    # Create an array of temperatures
    # Choose a step size for the temperature, dT = 0.1
    numTemps = int(3 / 0.1)
    temps = np.zeros(numTemps)
    temps[0] = 1
    dT = 0.1
    for i in range(1, numTemps):
        temps[i] = temps[i - 1] + dT

    # Create an arary's for the exact data
    exactMagArray = np.zeros(numTemps)
    exactHeatCapacity = np.zeros(numTemps)
    # Create array's for the exp data
    avgMagArray = np.zeros(numTemps)
    magneticSusceptibilityArray = np.zeros(numTemps)
    heatCapacityArray = np.zeros(numTemps)
    # Create array's for the uncertainties
    uncertMagArrray = np.zeros(numTemps)
    uncertMagSusArray = np.zeros(numTemps)
    uncertHeatCapArray = np.zeros(numTemps)
    # For a series of system sizes and temperatures run the simulation
    sizes = [8]
    for L in sizes:
        for i in range(0, numTemps):
            print("@ temp: ", temps[i])
            print("@ L = : ", L)
            # Experimental data
            avgMagArray[i] = temp.calculateAvgMagnetization()
            magneticSusceptibilityArray[i] = temp.calculateMagneticSusceptibilityPerSite()
            heatCapacityArray[i] = temp.calculateHeatCapcaityPerSite()
            # Exact data
            exactMagArray[i] = temp.calculateExactMagneticSusceptibilityPerSite()
            #print("avg energy: ", exactMagArray[i])
            exactHeatCapacity[i] = temp.calculateExactHeatCapacity()
            #print("heat capacity: ", exactHeatCapacity[i])
            # Uncertainties
            uncertMagArrray[i] = temp.bootStrappingMethod(temp.avgMagArray)
            uncertMagSusArray[i] = temp.bootStrappingMethod(magneticSusceptibilityArray)
            uncertHeatCapArray[i] = temp.bootStrappingMethod(temp.energyArray)
            print("The correlation time is: ", temp.correlationTime)
            print("The last three avg's were: ", temp.avgMagArray[1000], " ", temp.avgMagArray[1500], " ", temp.avgMagArray[2499], " ")
            #print("temps: ", temps)
        print("avgMagArray: ", avgMagArray)
        print("magneticSusceptibilityArray: ", magneticSusceptibilityArray)
        print("heatCapacityArray: ", heatCapacityArray)
        print("uncertaintyArrary: ", uncertMagArrray)
        print("uncertMagSusArray: ", uncertMagSusArray)
        print("uncertHeatCapArray: ", uncertHeatCapArray)
        plt.errorbar(temps, avgMagArray, uncertMagArrray)
        plt.plot(temps, exactMagArray)
        plt.title("Average Magnetization Per Site vs. Temperature")
        plt.xlabel("T")
        plt.ylabel("<|m|>")
        fileName = str(L) + "x" + str(L) + "AvgMag.pdf"
        plt.savefig(fileName)
        plt.close()
        plt.errorbar(temps, magneticSusceptibilityArray, uncertMagSusArray)
        plt.title("Magnetic Susceptibility Per Site vs. Temperature")
        plt.xlabel("T")
        plt.ylabel("x")
        fileName = str(L) + "x" + str(L) + "MagSuscept.pdf"
        plt.savefig(fileName)
        plt.close()
        plt.errorbar(temps, heatCapacityArray, uncertHeatCapArray)
        plt.plot(temps, exactHeatCapacity)
        plt.title("Heat Capacity Per Site vs. Temperature")
        plt.xlabel("T")
        plt.ylabel("Cv")
        fileName = str(L) + "x" + str(L) + "HeatCapacity.pdf"
        plt.savefig(fileName)
        plt.close()






#################
######TEST#######
#################
#mD4T1 = MonteCarlo(4,1)
#mD8T1 = MonteCarlo(8,1)
#mD16T1 = MonteCarlo(16,1)
#mD32T1 = MonteCarlo(32,1)

#mDT1 = MonteCarlo(4, 1)
#mDT1.metropolisAlgorithm(307200)
#mDT2 = MonteCarlo(8, 1)
#mDT2.metropolisAlgorithm(307200)
#mDT3 = MonteCarlo(16, 1)
#mDT3.metropolisAlgorithm(307200)
#mDT4 = MonteCarlo(32, 1)
#mDT4.metropolisAlgorithm(3000000)
simulation()
#mDT1.autocorr(mDT1.avgMagArray)
#mDT1.estimated_autocorrelation(mDT1.avgMagArray)
#mDT10 = MonteCarlo(32, 2.3)
#mDT10.metropolisAlgorithm(507200)

