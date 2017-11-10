# Xavier Linn
# MSE 215, Project 3
# A Monte Carlo Problem
# Metropolis algorithm to compute the properties of a 2D Ising models.

# Import statements
import numpy as np
import matplotlib.pyplot as plt

# Generate random elements
# Create a periodic indexing function
# etc..
# Create a grid class
class TwoDArray:
    # Notes:
    # array[row][col] = row[i][j]
    # Class variables out here. Shared by all TwoDArray instances.

    # Instantiate a 2D array
    def __init__(self, size):
        # Instance variables. Unique to each instance.
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
                if (randomValue < 6):
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
        #print(site)
        i = site[0]
        j = site[1]
        currentState = self.alloy[i][j]
        if currentState == 1:
            self.alloy[i][j] = -1
        else:
            self.alloy[i][j] = 1

# This class will run the simulation
class MonteCarlo:

    # Instantiate a  Monte Carlo Simulation
    def __init__(self, s, temp):
        self.twoDArray = TwoDArray(s)
        self.totalEnergy = 0
        self.size = s
        self.time = 0
        self.temperature = temp
        self.boltzmanConstant = 1
        self.beta = 1 / (self.temperature * self.boltzmanConstant)
        self.avgMagArray = np.array
        self.energyArray = np.array
        self.steps = 0
        self.ntrials = 0
        self.correlationTime = 0
        self.correlationFxn = np.array

    # Plot total energy vs time
    def totalEnergyVsTime(self):
        #print("Should be plotting now")
        plt.plot(self.energyArray)
        plt.show()

    # Plot avg mag vs time
    def magVsTime(self):
        #print("Should be plotting now")
        plt.plot(self.avgMagArray)
        plt.show()

    # Calculate the initial total energy, i.e. before any spins have been flipped
    def calculateInitialTotalEnergy(self):
        for i in range(0, self.size):
            for j in range(0, self.size):
                self.totalEnergy += self.twoDArray.alloy[i][j]*(self.twoDArray.stateAtIndexI(i + 1, j) + self.twoDArray.stateAtIndexJ(i, j + 1))
        self.totalEnergy *= (-1)

    # Calculate the initial total energy, i.e. before any spins have been flipped
    def calculateTotalEnergy(self):
        energy = 0
        for i in range(0, self.size):
            for j in range(0, self.size):
                energy += self.twoDArray.alloy[i][j]*(self.twoDArray.stateAtIndexI(i + 1, j) + self.twoDArray.stateAtIndexJ(i, j + 1))
        return energy * (-1)

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
        #self.twoDArray.printArray()
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
        sum = abs((np.sum(self.twoDArray.alloy, None, float) + 0.0) / np.power(self.size, 2))
        return sum

    # Metropolis Alogrithm
    def metropolisAlgorithm(self, steps):
        print("Running MD Simulation")
        print("size: ", self.size)
        print("steps: ", steps)

        # Setup data structures
        numSpin = self.size * self.size # By convention one "time" step is L^2 spin flip attempts
        initialNumDataPoints = int(steps / numSpin)
        #print("initialNumDataPoints: ", initialNumDataPoints)
        self.avgMagArray = np.arange(initialNumDataPoints, dtype=np.float)
        self.energyArray = np.arange(initialNumDataPoints)

        # 0. Take the system to equilibrium.

        self.calculateInitialTotalEnergy() # Calculate the total energy before any spin flip attempts
        for i in range(0, steps):
            site = self.generateRandomSite() # Select a random site
            totalEnergyAfterSpinFlip = self.calculateTotalEnergyAfterSpinFlip(site) # Calculate the energy after a random spin is flipped
            # 3. Accept or reject the spin flip.
            flip = self.acceptFlipOrNot(totalEnergyAfterSpinFlip)
            if (flip):
                self.totalEnergy = totalEnergyAfterSpinFlip
            else:
                self.twoDArray.changeState(site)
            # Store energy and magnetization while system reaches equilibrium
            if (i % numSpin == 0):
                if (i / numSpin == initialNumDataPoints):
                    #print(i)
                    #print("Size of avgMagArray w/ non-eq. data: ", self.avgMagArray.size)
                    #print("Size of energyArray w/ non-eq. data: ", self.energyArray.size)
                    #print("self.time: ", self.time)
                    break
                self.avgMagArray[self.time] = self.calculateAvgMagnetization()
                self.energyArray[self.time] = self.totalEnergy
                self.time += 1

        # Plot before removing non-equilibrium data
        #self.totalEnergyVsTime()
        #self.magVsTime()


        # Determine the correlation time.
        self.correlationFxn = self.estimated_autocorrelation(self.avgMagArray)
        self.correlationTime = self.integrateCorrelation(self.correlationFxn)
        #print("Size of avgMagArray w/ non-eq. data: ", self.avgMagArray.size)
        #print("Size of energyArray w/ non-eq. data: ", self.energyArray.size)
        # Remove non-equilibrium data using correlation time. "Collect Data" now that the system is at equilibrium
        scaledCorrelationTime = int(2 * self.correlationTime)
        #print("Correlation time; ", self.correlationTime)
        #print("Correlation time scaled; ", scaledCorrelationTime)
        self.avgMagArray = self.avgMagArray[scaledCorrelationTime:]
        self.energyArray = self.energyArray[scaledCorrelationTime:]
        #print("Size of avgMagArray w/ only eq. data: ", self.avgMagArray.size)
        #print("Size of energyArray w/ only eq. data: ", self.energyArray.size)

        # Determine the uncertainty
        #uncertInMag = self.bootStrappingMethod()
        #print("The uncertainty in avg mag is: ", uncertInMag)

        # Plot after removing non-equilibrium data
        #self.totalEnergyVsTime()
        #self.magVsTime()

    def estimated_autocorrelation(self, x):
        n = len(x)
        variance = x.var()
        x = x-x.mean()
        r = np.correlate(x, x, mode = 'full')[-n:]
        #assert N.allclose(r, N.array([(x[:n-k]*x[-(n-k):]).sum() for k in range(n)]))
        result = r/(variance*(np.arange(n, 0, -1)))
        #print("THE ZERO IS APPROX AT INDEX ", self.findZeroIndex(result))
        #plt.plot(result)
        #plt.show()
        return result

    def findZeroIndex(self, x):
        i = 0
        j = 1

        #if (self.size == 4):
        #    return 2

        for v in x:
            #print("v: ", v)
            #print("x[i]: ", x[i])
            #print("x[j]: ", x[j])
            if x[i] > 0 and x[j] < 0:
                #print("THE VALUE IS: ", v)
                #print("THE INDEX IS: ", i)
                return j
            i += 1
            j += 1

    def integrateCorrelation(self, correlationFxn):
        index = self.findZeroIndex(correlationFxn)
        correlationFxn = correlationFxn[:index + 1]
        #plt.plot(correlationFxn)
        #plt.show()
        correlationTime = np.trapz(correlationFxn)
        return correlationTime

    def calculateNumberOfIndpTrials(self):
        #print(self.avgMagArray.size)
        self.ntrials = int(self.avgMagArray.size / (2 * self.correlationTime))

    # Calculates the the variance 
    def calculateHeatCapcaityPerSite(self):
        print("Calculating head capacity")
        # Cv = (1 / (Kb * T^2 * L^2)) * (<E^2> - <E>^2)
        var = np.var(self.energyArray)
        Cv = (1 / (np.power(self.size, 2) * np.power(self.temperature, 2))) * var


    # Caclulate the magnetic susceptibility per site 
    def calculateMagneticSusceptibilityPerSite(self):
        print("Calculating Magnetic Susceptibility")
        # x = beta * (1 / L^2) (<m^2> - <m>^2), beta = 1 / T
        var = np.var(self.avgMagArray)
        x = self.beta * (1 / np.power(self.size, 2)) * var



    def bootStrappingMethod(self, data):
        self.calculateNumberOfIndpTrials() # Calculate the number in indep. trials
        #print("ntrials: ", self.ntrials)
        p = np.arange(self.ntrials, dtype='f')
        #print("initial p is: ", p)
        #print("ntrials: ", self.ntrials)
        #n = int(np.power(self.ntrials, 2))
        n = 1000
        avgSampleArray = np.arange(n, dtype='f')
        #print(n)
        for i in range(0, n):
            for j in range(0, int(self.ntrials)):
                index = np.random.randint(1, data.size)
                #print("index: ", index)
                #print("avgmag ", self.avgMagArray[index])
                #print("j: ", j)
                p[j] = data[index]
                #p[j] = self.avgMagArray[index]
            #print(p)
            avgSample = p.mean()
            #print("avgSample: ", avgSample)
            avgSampleArray[i] = avgSample
        #print(avgSampleArray)
        avp = np.sqrt(np.var(avgSampleArray))
        return avp

def calculateMagVsTemp():
    # Create an array of temperatures
    # Choose a step size for the temperature, try dT = 0.1
    numTemps = int(3 / 0.1)
    print("numTemps: ", numTemps)
    temps = np.zeros(numTemps)
    print("temps: ", temps)
    #temps.fill(0)
    temps[0] = 1
    dT = 0.1
    for i in range(1, numTemps):
        temps[i] = temps[i - 1] + dT
    print(temps)

    # Create an array of error bars for each temperature
    # Create an array for the avg mag at each temperature
    avgMagArray = np.zeros(numTemps)
    magneticSusceptibilityArray = np.zeros(numTemps)
    heatCapacityArray = np.zeros(numTemps)
    #avgMagArray.fill(0)
    # Create an array for the uncertainty
    uncertMagArrray = np.zeros(numTemps)
    uncertMagSusArray = np.zeros(numTemps)
    uncertHeatCapArray = np.zeros(numTemps)
    #uncertaintyArrray.fill(0)
    sizes = [4, 8, 16, 32]
    for L in sizes:
        for i in range(0, numTemps):
            print("@ temp: ", temps[i])
            print("@ L = : ", L)
            temp = MonteCarlo(L, temps[i])
            temp.metropolisAlgorithm(3000000)
            avgMagArray[i] = temp.calculateAvgMagnetization()
            magneticSusceptibilityArray[i] = temp.calculateMagneticSusceptibilityPerSite()
            heatCapacityArray[i] = temp.calculateHeatCapcaityPerSite()
            uncertMagArrray[i] = temp.bootStrappingMethod(avgMagArray)
            uncertMagSusArray[i] = temp.bootStrappingMethod(magneticSusceptibilityArray)
            uncertHeatCapArray[i] = temp.bootStrappingMethod(heatCapacityArray)
            plt.errorbar(temps, avgMagArray, uncertMagArrray)
            plt.show()
            plt.errorbar(temps, magneticSusceptibilityArray, uncertMagSusArray)
            plt.show()
            plt.errorbar(temps, heatCapacityArray, uncertHeatCapArray)
            plt.show()  
            plt.show()
        print("temps: ", temps)
        print("avgMagArray: ", avgMagArray)
        print("uncertaintyArrary: ", uncertMagArrray)
    plt.errorbar(temps, avgMagArray, uncertMagArrray)
    plt.show()
    plt.errorbar(temps, magneticSusceptibilityArray, uncertMagSusArray)
    plt.show()
    plt.errorbar(temps, heatCapacityArray, uncertHeatCapArray)
    plt.show()  
    plt.show()




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
calculateMagVsTemp()
#mDT1.autocorr(mDT1.avgMagArray)
#mDT1.estimated_autocorrelation(mDT1.avgMagArray)
#mDT10 = MonteCarlo(32, 10)
#mDT10.metropolisAlgorithm(307200)

