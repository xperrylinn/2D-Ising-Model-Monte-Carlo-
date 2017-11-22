# Xavier Linn
# MSE215, Fall 2017
# LAMMPS Molecualr Dynamics
import re
import numpy as np
import matplotlib.pyplot as plt

# ##########PRACTICE############
# # m = re.match(r"(\w+) (\w+)", "Isaac Newton, physicist")
# m = re.match("(?:\s*)(-?[0-9]*\.?[0-9]*)(?:\s+)(-?[0-9]*\.?[0-9]*)(?:\s+)(-?[0-9]*\.?[0-9]*)(?:\s+)(-?[0-9]*\.?[0-9]*)(?:\s+)(-?[0-9]*\.?[0-9]*)(?:\s*)", " 0         1112   -1874.2011    -1945.926    71.724903, physicist")
# print(m.group(0), "0")
# print(m.group(1), "1")
# print(m.group(2), "2")
# print(m.group(3), "3")
# print(m.group(4), "4")
# print(m.group(5), "5")

##############################

############PART1################
###########QUESTION1#############
print("Hello World: Question 1")

# Process the LAMMPS std.out file
stdOutFile = open("stdoutPd.txt", "r") # Open the file
lines = stdOutFile.readlines() # Put all lines of the file into a list
lines = lines[8:] # Remove the first 8 lines, they are part of the header
lines = lines[:-19] # Remove the last 19 lines, not needed

# Create arrays to hold desired data
numberDataPoints = len(lines)
steps = np.arange(numberDataPoints)
temp = np.arange(numberDataPoints)
totalEnergy = np.arange(numberDataPoints)
potentialEnergy = np.arange(numberDataPoints)
kineticEnergy = np.arange(numberDataPoints)

# Store the LAMMPS data
i = 0 # index, used to set data of data array's
for line in lines:
    # Order of data in file: Step Temp TotalEnergy PotentialEnergy KineticEnergy
    m = re.match("(?:\s*)(-?[0-9]*\.?[0-9]*)(?:\s+)(-?[0-9]*\.?[0-9]*)(?:\s+)(-?[0-9]*\.?[0-9]*)(?:\s+)(-?[0-9]*\.?[0-9]*)(?:\s+)(-?[0-9]*\.?[0-9]*)(?:\s*)", line) # Python Regex, parse the input file for desired info
    steps[i] = float(m.group(1))
    temp[i] = float(m.group(2))
    totalEnergy[i] = float(m.group(3))
    potentialEnergy[i] = float(m.group(4))
    kineticEnergy[i] = float(m.group(5))
    i += 1

# # Plot total energy vs. steps
# plt.title("Total Energy (eV) vs. Steps")
# plt.xlabel("Steps")
# plt.ylabel("Total Energy")
# plt.grid('on')
# plt.plot(steps, totalEnergy)
# plt.show()
# plt.close()
# # Plot kinetic energy vs. steps
# plt.title("Potential Energy (eV) vs. Steps")
# plt.xlabel("Steps")
# plt.ylabel("Potential Energy")
# plt.grid('on')
# plt.plot(steps, potentialEnergy)
# plt.show()
# plt.close()
# # Plot kinetic energy vs. steps
# plt.title("Kinetic Energy (eV) vs. Steps")
# plt.xlabel("Steps")
# plt.ylabel("Kinetic Energy")
# plt.grid('on')
# plt.plot(steps, kineticEnergy)
# plt.show()
# plt.close()

############PART2###############
###########QUESTION2#############
print("Hello World: Question 2")

# TODO: process velocities dump file

##########PRACTICE############
# m = re.match(r"(\w+) (\w+)", "Isaac Newton, physicist")
m2 = re.match("(?:\s*)(-?[0-9]*\.?[0-9]*)(?:\s+)(-?[0-9]*\.?[0-9]*)(?:\s+)(-?[0-9]*\.?[0-9]*)(?:\s+)(-?[0-9]*\.?[0-9]*)(?:\s+)(-?[0-9]*\.?[0-9]*)(?:\s*)", "ITEM: ATOMS id type vx vy vz")
if (m2 is None):
    print("NONE - NO MATCH")
else:
    print(m2.group(0), "ground(0)")
    print(m2.group(1), "group(1)")
    print(m2.group(2), "group(2)")
    print(m2.group(3), "group(3)")
    print(m2.group(4), "group(4)")
    print(m2.group(5), "group(5)")

# Process the LAMMPS velocities.dump file
velocitiesFile = open("velocitiesdump2.txt", "r") # Open the file
lines2 = velocitiesFile.readlines() # Put all lines of the file into a list

#Test
print("The size of lines2 is: ", len(lines2))
print("MATH: ", len(lines2) - 9 * 100)

# Create arrays to hold the data
numberDataPoints2 = 500 * 100 # 500 velocities per coordinate, have 100 measurements
# atomNumber = np.arange(numberDataPoints2)
# velocities = np.arange(numberDataPoints2)
velocities = []

# Store the LAMMPS data
i = 0 # index, used to set data of data array's
for line in lines2:
    # Order of data: id type vx vy vz
    m2 = re.match("(?:\s*)(-?[0-9]*\.?[0-9]*)(?:\s+)(-?[0-9]*\.?[0-9]*)(?:\s+)(-?[0-9]*\.?[0-9]*)(?:\s+)(-?[0-9]*\.?[0-9]*)(?:\s+)(-?[0-9]*\.?[0-9]*)(?:\s*)", line)
    if (m2 is None): #
        # Do nothing
        a = None
    else:
        # velocities[i] = float(m2.group(3)) # Take only group 3, the vx
        velocities.append(float(m2.group(3))) # Take only group 3, the vx

    i += 1

print("velocities has size: ", len(velocities))
print(velocities)



