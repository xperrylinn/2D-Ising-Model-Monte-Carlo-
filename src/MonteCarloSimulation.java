/**
 * Created by xperrylinn on 10/17/17.
 */

// Imports
import sun.jvmstat.monitor.MonitorException;

import java.util.Random;

/*
    TODO: Create a separate classes. One to run the simulation, one to do the
    math, one to create the grid
 */
public class MonteCarloSimulation {
    int[][] alloy;
    int size;
    Random rand;

    /*
        Class constructor.
     */
    public MonteCarloSimulation(int s) {
        size = s;   // square grid, with length size
        alloy = new int[size][size]; // 2D array data structure
        rand = new Random();    // instantiate a random number generator
    }

    /*
        Populate the 2D array randomly with 1's and 0's
     */
    public void populateArray() {
        for (int i = 0; i < size; i += 1) {
            for (int j = 0; j < size; j += 1) {
                double r = rand.nextDouble(); // Generate a random number, r, between 0.0 and 1.0
                if (r > 0.5) { // if r it less than 0.5, assign a 1 to the grid point at i, j
                    alloy[i][j] = 1;
                } else { // if r it less than 0.5, assign a 0 to the grid point at i, j
                    alloy[i][j] = 0;
                }
            }
        }
    }

    /*
        Print the 2D array
     */
    public void printArray() {
        System.out.println("Printing the 2D array");
        for (int i = 0; i < size; i += 1) {
            for (int j = 0; j < size; j += 1) {
                System.out.print(alloy[i][j] + " ");
            }
            System.out.println();
        }
    }

    /*
        TODO: Create an index function, periodic boundry conditions
        from 0 to L - 1. Given L, return 0. Given -1 give L - 1.
     */
    public void getElement(int i, int j) {
        if (false) { // If i, j is in the boundary size x size

        } else if (false) {

        } else {

        }
    }

    public static void main(String[] args) {
        // Prints "Hello, World" to the terminal window.
        System.out.println("MonteCarloSimulation");
        // Create a new MonteCarloSimulation
        MonteCarloSimulation test = new MonteCarloSimulation(5);
        test.populateArray();
        test.printArray();
    }
}
