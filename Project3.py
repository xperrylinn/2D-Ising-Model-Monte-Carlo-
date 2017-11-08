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
        #fileE = open("energyVsTime.txt", "w")
        self.avgMagArray = np.arange(steps / 300, dtype=np.float)
        #fileM = open("avgMagVsTime.txt", "w")
        self.energyArray = np.arange(steps / 300)
        #self.twoDArray.printArray()
        # 1. Establish an initial configuration of moments.
        self.calculateInitialTotalEnergy()
        for i in range(0, steps):
            # 2. Select a random moment. Compute change in change
            #    in energy associated with flipping the spin.
            site = self.generateRandomSite()
            #self.twoDArray.changeState(site)
            #newTotalEnergyEnergy = self.calculateTotalEnergy()
            newTotalEnergyEnergy = self.calculateTotalEnergyAfterSpinFlip(site)
            # 3. Accept or reject the spin flip.
            flip = self.acceptFlipOrNot(newTotalEnergyEnergy)
            if (flip):
                self.totalEnergy = newTotalEnergyEnergy
            else:
                self.twoDArray.changeState(site)
            if (i % 300 == 0):
                self.avgMagArray[self.time] = self.calculateAvgMagnetization()
                self.energyArray[self.time] = self.totalEnergy
                self.time += 1
            # 4. Compute quantities of interest: Magnetization, Total Energy
            # 5. Repeat from step 2 as needed.
        # 6. Analyze the results of the simulation.
        self.correlationFxn = self.estimated_autocorrelation(self.avgMagArray)
        self.correlationTime = self.integrateCorrelation(self.correlationFxn)
        print("The correlation time is: ", self.correlationTime)
        print("The uncertainty in avg mag is: ", self.bootStrappingMethod())
        #self.totalEnergyVsTime()
        #self.magVsTime()

    def autocorr(self, x):
        result = np.correlate(x, x, 'full')
        #plt.plot(result)
        #plt.show()
        #return result[result.size/2:]

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
        self.ntrials = self.time / (2 * self.correlationTime)

    def bootStrappingMethod(self):
        self.calculateNumberOfIndpTrials()
        p = np.arange(int(self.ntrials), dtype='f')
        #print("initial p is: ", p)
        #print("ntrials: ", self.ntrials)
        #n = int(np.power(self.ntrials, 2))
        n = 1000
        avgSampleArray = np.arange(n, dtype='f')
        #print(n)
        for i in range(0, n):
            for j in range(0, int(self.ntrials)):
                index = np.random.randint(1, self.time)
                #print("index: ", index)
                #print("avgmag ", self.avgMagArray[index])
                #print("j: ", j)
                p[j] = self.avgMagArray[index]
            #print(p)
            avgSample = p.mean()
            #print("avgSample: ", avgSample)
            avgSampleArray[i] = avgSample
        #print(avgSampleArray)
        avp = np.sqrt(np.var(avgSampleArray))
        return avp



#################
######TEST#######
#################
#mD4T1 = MonteCarlo(4,1)
#mD8T1 = MonteCarlo(8,1)
#mD16T1 = MonteCarlo(16,1)
#mD32T1 = MonteCarlo(32,1)

mDT1 = MonteCarlo(4, 1)
mDT1.metropolisAlgorithm(307200)
mDT2 = MonteCarlo(8, 1)
mDT2.metropolisAlgorithm(307200)
mDT3 = MonteCarlo(16, 1)
mDT3.metropolisAlgorithm(307200)
mDT4 = MonteCarlo(32, 1)
mDT4.metropolisAlgorithm(307200)
#mDT1.autocorr(mDT1.avgMagArray)
#mDT1.estimated_autocorrelation(mDT1.avgMagArray)
#mDT10 = MonteCarlo(32, 10)
#mDT10.metropolisAlgorithm(307200)

