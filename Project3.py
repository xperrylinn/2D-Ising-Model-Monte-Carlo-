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
class TwoDArrary:
    # Notes:
    # array[row][col] = row[i][j]
    # Class variables out here. Shared by all TwoDArray instances.

    #
    def __init__(self, size):
        # Instance variables. Unique to each instance.
        self.size = size # 2D array is a square with length size. 0, 1, ... size - 1
        self.alloy = np.random.randint(0,2,(size,size)) # Create a 2D array with populated randomly with 1's and 0's

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

    # Retrieve the state at index i
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
            self.alloy[i][j] = 0
        else:
            self.alloy[i][j] = 1

# Create a test TwoDArray
test = TwoDArrary(4)
test.printArray()
#print("Testing indexing i, on edge cases")
#print("(-1, 0):", test.stateAtIndexI(-1,0))
#print("(4, 0):", test.stateAtIndexI(4,0))
#print("Testing indexing j, on edge cases")
#print("(0, -1):", test.stateAtIndexJ(0,-1))
#print("(0, 4):", test.stateAtIndexJ(0,4))
print("Now change the state at index i, j = 1, 1")
print("before: ", test.alloy[1][1])
print("Calling change state..")
test.changeState(1,1)
print("after: ", test.alloy[1][1])
test.printArray()
