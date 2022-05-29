"""
Velocity Verlet Time Integrator to simulate the evolution of a solar system of Celestial bodies, described in the Body3D class.

where the parameters are read in from a data file
and passed to the functions that
calculate force and potential energy.

Using units that make for better understanding of the simulated system.

G = 8.8877 × 10−10 AU^3 M_earth ^ -1 day ^ -2

Author : E Forster, s1639706
"""

import sys
import math
import numpy as np
import matplotlib.pyplot as pyplot
from scipy.signal import find_peaks
from body3D import Body3D

def calculate_body_separation(body_list) :

    N = len(body_list)
    separations_matrix = np.zeros((N, N, 3))

    for body in body_list :
        for i in range(N) :
            for j in range(i + 1, N) :

                separation = body_list[i].position - body_list[j].position

                separations_matrix[i, j] = separation
                separations_matrix[j, i] = - separation  # Applies symmetry to reduce computing time

    return separations_matrix


def calculate_gravitational_force(body_list, G) :

    N = len(body_list)
    separation_matrix = calculate_body_separation(body_list)

    force_matrix = np.zeros((N, N, 3))

    for i in range(N):
        for j in range(i + 1, N):

            modulus_sep_matrix = np.linalg.norm(separation_matrix[i, j])

            if modulus_sep_matrix == 0:  # Applies cut off radius to reduce computing time

                force_matrix[i, j] = 0

            else:
                force_matrix[i, j] = - G * body_list[i].mass * body_list[j].mass * (separation_matrix[i, j] / (modulus_sep_matrix ** 3))
                force_matrix[j, i] = - force_matrix[i, j]  # Applies Newton's 3rd Law to reduce computing time

    return force_matrix

def calculate_gravitational_potential(body_list, G) :

    N = len(body_list)
    separation_matrix = calculate_body_separation(body_list)
    g_potential = 0

    for i in range(N) :
        for j in range(i + 1, N) :

            modulus_separation_matrix = np.linalg.norm(separation_matrix[i, j])

            g_potential += - G * body_list[i].mass * body_list[j].mass / modulus_separation_matrix

    return g_potential

"""
def oscillations(pos_list, dt) :

    Method to calculate the period of the oscillating body

    
    wave_peaks = find_peaks(pos_list)                           # Uses scipy.signal find_peaks function to index all peaks
    peak_pos = wave_peaks[0]                                    # Creates an array of the indices of all peaks

    first_peak = peak_pos[0]                                    # First peak
    second_peak = peak_pos[1]                                   # Second peak

    period = (len(pos_list[first_peak : second_peak])) * dt     # Calculates period

    return period
"""

# Begin main code
def main() :
    """
    The main method carries out the simulation in a few parts:

    1.) Reads in data file from the command line
    2.) Specifies initial conditions
    3.) Initialises data lists for plotting later
    4.) Starts a time integration loop
    5.) Plots particle trajectory to screen
    6.) Plots particle energy to screen
    7.) Measures the frequency of oscillations and prints the wave-number and it's inaccuracy to screen
    8.) Measures the energy inaccuracy of the simulation and prints it to the screen

    """
    with open(sys.argv[1], "r") as infile :

        # Part 1.) Reads in data file from the command line

        # Read name of files from command line which needs 3 parts
        if len(sys.argv) != 2 :

            # Helpful error message if the format is incorrect
            print("Wrong number of arguments.")
            print("Usage: " + sys.argv[0] + "<input file>")
            quit()

        else :

            line1 = infile.readline()           # Processes line 1 of input file
            line2 = infile.readline()
            line2 = line2.split()

            if len(line2) != 2 :

                print("Wrong number of arguments in line 2 of input file, dt and numstep. ")

            else :

                dt = line2[0]
                numstep = line2[1]

            line3 =  infile.readline()
            line4 = infile.readline()
            line4 = line4.split()

            if len(line4) != 2 :

                print("Wrong number of arguments in line 4 of input file, trajectory outfile and energy outfile. ")

            else :

                trajectory_xyz = line4[0]
                energy_output = line4[1]

            line5 = infile.readline()
            line6 = infile.readline()
            line6 = line6.split()

            if len(line6) != 1 :

                print("Wrong number of arguments in line 6 of input file, Body data. ")

            else :

                body_data = line6[0]

            # Open output file
            trajectories = open(trajectory_xyz, "w")
            energies = open(energy_output, "w")

            energies.write("Time, Potential Energy, Kinetic Energy, Total Energy \n")

    infile.close()

    # Part 2.) Specifies initial conditions

    time = 0.0
    G = 8.8877 * 10 ** - 10
    body_list = []

    bodies = np.genfromtxt(body_data, dtype='str')
    N = len(bodies)

    for body in range(N) :

        body_list.append(Body3D(label = f"Body_{body}", mass = 1, position = np.zeros(3), velocity = np.zeros(3)))
        print(body_list)
        exit()

    # Get initial force
    force_matrix = calculate_gravitational_force(body_list, G)


    # Part 3.) Initialises data lists for plotting later

    time_list = [time]


    # Part 4.) Starts a time integration loop

    for i in range(numstep) :

        # Update body position


        # Update force


        # Update body velocity by averaging current and new forces
        
        # Re-define force value


        # Increase time
        time += dt
        
        # Output body information



        # Append information to data lists
        time_list.append(time)
        energy_list.append(energy)
    

    # Post-simulation:
    # Close output file
    outfile.close()

    # Part 5.) Plots body trajectory to screen

    # Part 6.) Plots body energy to screen

    # Plot body energy to screen
    pyplot.title('Velocity Verlet : Total Energy vs Time')
    pyplot.xlabel('Time : ')
    pyplot.ylabel('Energy : ')
    pyplot.plot(time_list, energy_list)
    pyplot.show()

    # Part 8.) Measures the energy inaccuracy of the simulation and prints it to the screen

    initial_energy = total_energy_list[0]
    max_energy = max(total_energy_list)
    min_energy = min(total_energy_list)

    delta_energy = max_energy - min_energy
    energy_inaccuracy = delta_energy / initial_energy

    print("Energy inaccuracy : +/-", energy_inaccuracy, "[energy] ")

# Execute main method, but only when directly invoked
if __name__ == "__main__" :
    main()
