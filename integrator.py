"""
Velocity Verlet Time Integrator to simulate the evolution of a solar system of Celestial bodies, described in the Body3D class.

where the parameters are read in from a data file
and passed to the functions that
calculate force and potential energy.

Author : E Forster, s1639706
"""

import sys
import math
import numpy as np
import matplotlib.pyplot as pyplot
from scipy.signal import find_peaks
from body3D import Body3D

# G = 8.8877 × 10−10 AU^3 M_earth ^ -1 day ^ -2

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
                force_matrix[i, j] = - G * body_list[i].mass * body_list[j].mass * (separation_matrix[i, j] / (np.mod(separation_matrix[i, j]) ** 3))
                force_matrix[j, i] = - force_matrix[i, j]  # Applies Newton's 3rd Law to reduce computing time

    return force_matrix



def oscillations(pos_list, dt) :
    """
    Method to calculate the period of the oscillating particle
    and then obtain the wave-number given by:

    wave-number = 1 / wavelength

    :param pos_list: position list along the time-axis
    :param dt: time-step for simulation

    :return: wave_number
    """
    wave_peaks = find_peaks(pos_list)                           # Uses scipy.signal find_peaks function to index all peaks
    peak_pos = wave_peaks[0]                                    # Creates an array of the indices of all peaks

    first_peak = peak_pos[0]                                    # First peak
    second_peak = peak_pos[1]                                   # Second peak

    period = (len(pos_list[first_peak : second_peak])) * dt     # Calculates period
    period_in_seconds = period * 1                # Converts period to seconds

    frequency = 1 / period_in_seconds                           # Units : Hz
    speed_of_light = 299792458                                  # Units : m/s
    wavelength = (speed_of_light / frequency) * 100             # Units : m (multiply by 100 to convert to cm) -> cm

    wave_number = 1 / wavelength
    return wave_number                                          # Units : cm^-1

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
        if len(sys.argv) != 3 :

            # Helpful error message if the format is incorrect
            print("Wrong number of arguments.")
            print("Usage: " + sys.argv[0] + "<input file>" + "<output file>")
            quit()
        else :
            outfile_name = sys.argv[2]

            # Open output file
            outfile = open(outfile_name, "w")

            line1 = infile.readline()           # Processes line 1 of input file


    infile.close()

    # Part 2.) Specifies initial conditions

    time = 0.0
    G = 8.8877 * 10 ** - 10

    # Get initial force
    force1 = []
    force2 = - force1
    energy = []

    # Part 3.) Initialises data lists for plotting later

    time_list = [time]
    position_list = []
    energy_list = [energy]

    # Part 4.) Starts a time integration loop

    for i in range(numstep) :

        # Update body position
        p1.update_2nd_position(dt, force1)

        # Update force
        force1_new = []
        force2_new = - force1_new

        # Update body velocity by averaging current and new forces
        
        # Re-define force value
        force1 = force1_new
        force2 = force2_new

        # Increase time
        time += dt
        
        # Output body information
        energy = p1.kinetic_e()


        # Append information to data lists
        time_list.append(time)
        energy_list.append(energy)
    

    # Post-simulation:
    # Close output file
    outfile.close()

    # Part 5.) Plots body trajectory to screen

    pyplot.title('Velocity Verlet : Position vs Time')
    pyplot.xlabel('Time : ')
    pyplot.ylabel('Position : ')
    pyplot.plot(time_list)
    pyplot.plot(time_list )
    pyplot.show()

    # Part 6.) Plots body energy to screen

    # Plot body energy to screen
    pyplot.title('Velocity Verlet : Total Energy vs Time')
    pyplot.xlabel('Time : ')
    pyplot.ylabel('Energy : ')
    pyplot.plot(time_list, energy_list)
    pyplot.show()

    # Part 7.) Measures the frequency of oscillations and prints the wave-number and it's inaccuracy to screen

    wave_number = oscillations(position_list, dt)
    print("Wave-number :", wave_number, "cm ^ -1 .")
    v_nought = 1525.5847124925656
    delta_v = wave_number - v_nought
    wave_number_inaccuracy = delta_v / v_nought
    print("Wave-number inaccuracy : +/-", wave_number_inaccuracy, "cm ^ -1 .")

    # Part 8.) Measures the energy inaccuracy of the simulation and prints it to the screen

    initial_energy = energy_list[0]
    max_energy = max(energy_list)
    min_energy = min(energy_list)

    delta_energy = max_energy - min_energy
    energy_inaccuracy = delta_energy / initial_energy

    print("Energy inaccuracy : +/-", energy_inaccuracy, "eV ")


# Execute main method, but only when directly invoked
if __name__ == "__main__" :
    main()
