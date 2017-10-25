# Xavier Linn
# MSE 215, Project 3
# A Monte Carlo Problem
# Metropolis algorithm to compute the properties of a 2D Ising models.

# Import statements
import numpy as np

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
    def changeState(self, i, j):
        currentState = self.alloy[i][j]
        if currentState == 1:
            self.alloy[i][j] = -1
        else:
            self.alloy[i][j] = 1

# This class will run the simulation
class MonteCarloSimulation:

    # Instantiate a  Monte Carlo Simulation
    def __init__(self, s):
        self.twoDArray = TwoDArray(s)
        self.totalEnergy = 0
        self.size = s



    # Calculate the total energy
    # Toatl energt eq.
    # J*Sum(i = 0 to L-1)Sum(j = 0 to L -1)(s(i,j)s(i,j+1)+s(i,j)s(i+1,j)) - H*Sum(i = 0 to L-1)Sum(j = 0 to L -1)(s(i,j))
    # Set J = 1, H = 0
    def calculateTotalEnergy(self):
        print("Caclulating the total energy")
        # Iterate through the 2D array calculate
        # Sum up the energy in the totalEnergy cache variable
        #print(self.twoDArray.printArray())
        #print(self.totalEnergy)
        for i in range(0, self.size - 1):
            for j in range(0, self.size - 1):
                self.totalEnergy += self.twoDArray.alloy[i][j]*(self.twoDArray.stateAtIndexI(i + 1, j) + self.twoDArray.stateAtIndexJ(i, j + 1))
        #print(self.totalEnergy)


    # Calculate change in energy:
    # calculates the energy between two iterations
    def flipSpinCalcNewEnergy(self):
        print("Caclulating the change in energy")

        # Select a random site i, j
        i = np.random.randint(0, self.size)
        j = np.random.randint(0, self.size)

        #
        print("site: ",i ," ",j)

        #
        print("Total energy before flip", self.totalEnergy)

        #Calcualte energy contributution at side i,j before flip
        energyContributionBeforeFlip = self.twoDArray.alloy[i][j] * (self.twoDArray.stateAtIndexI(i + 1, j) + self.twoDArray.stateAtIndexJ(i, j + 1) + self.twoDArray.stateAtIndexI(i - 1, j) + self.twoDArray.stateAtIndexJ(i, j - 1))
        print("Energy contribution at i,j before flip: ", energyContributionBeforeFlip)

        # Subtract energy contribution at site i,j before flip from total energy
        self.totalEnergy -= energyContributionBeforeFlip

        #
        print("Printing array before splin flip")
        self.twoDArray.printArray()

        # Change the spin at a random site
        self.twoDArray.changeState(i, j)

        #
        print("Printing array after splin flip")
        self.twoDArray.printArray()

        energyContributionAfterFlip = self.twoDArray.alloy[i][j] * (self.twoDArray.stateAtIndexI(i + 1, j) + self.twoDArray.stateAtIndexJ(i, j + 1) + self.twoDArray.stateAtIndexI(i - 1, j) + self.twoDArray.stateAtIndexJ(i, j - 1))
        print("Energy contribution at i,j after flip: ", energyContributionAfterFlip)

        # Subtract energy contribution at site i,j before flip from total energy
        self.totalEnergy += energyContributionAfterFlip

        print("Total energy after flip", self.totalEnergy)



    # Accept the spin flip or not
    def acceptFlipOrNot(self):
        print("Determining to accept spin flip or not")

# Notes: Set J = 1 in equation 7
# Set H = 0 in equation 7
#
#




#################
######TEST#######
#################

# Create a test TwoDArray with random
#test = TwoDArrary(4)
#test.printArray()

# Test a totalEnergy calculation
mD = MonteCarloSimulation(6)
mD.calculateTotalEnergy()
mD.flipSpinCalcNewEnergy()

#print("Testing indexing i, on edge cases")
#print("(-1, 0):", test.stateAtIndexI(-1,0))
#print("(4, 0):", test.stateAtIndexI(4,0))
#print("Testing indexing j, on edge cases")
#print("(0, -1):", test.stateAtIndexJ(0,-1))
#print("(0, 4):", test.stateAtIndexJ(0,4))
#print("Now change the state at index i, j = 1, 1")
#print("before: ", test.alloy[1][1])
#print("Calling change state..")
#test.changeState(1,1)
#print("after: ", test.alloy[1][1])
#test.printArray()
