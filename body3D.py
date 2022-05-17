"""
 Body3D is a class describing a Celestial Body.

 An instance describes a Celestial Body in Euclidean 3D space:
 velocity and position are [3] arrays

 Includes time integrator methods + calculation methods + update methods

 Author: E Forster, s1639706

"""
import math
import numpy as np


class Body3D(object) :
    """
    Class to describe celestial bodies in 3D space

        Properties:
    label: name of the body
    mass: mass of the body
    pos: position of the body
    vel: velocity of the body

        Methods:
    __init__ - initialises a body in 3D space
    __str__ - sets up xyz of body
    kinetic_e  - computes the kinetic energy
    momentum - computes the linear momentum
    update_pos - updates the position to 1st order
    update_pos_2nd - updates the position to 2nd order
    update_vel - updates the velocity

        Static Methods:
    new_p3d - initializes a Body3D instance from a file handle
    sys_kinetic - computes total K.E. of a Body3D list
    com_vel - computes total mass and CoM velocity of a Body3D list
    """

    def __init__(self, label, mass, position, velocity) :
        """
        Initialises a celestial body in 3D space

        :param label: String w/ the name of the body
        :param mass: float, mass of the body
        :param position: [3] float array w/ position
        :param velocity: [3] float array w/ velocity
        """
        self.label = label
        self.mass = mass
        self.position = position
        self.velocity = velocity


    def __str__(self) :
        """
        XYZ-compliant string. The format is
        <label>    <x>  <y>  <z>

        :return xyz_string: (label, x, y, z)
        """
        xyz_string = str(self.label + "    " + str(self.position[0]) + "   " + str(self.position[1]) + "   " + str(self.position[2]))
        return xyz_string


    def kinetic_energy(self) :
        """
        Returns the kinetic energy of a Body3D instance

        :return ke: float, 1/2 m v**2
        """
        kinetic_energy = (1/2) * self.mass * (np.linalg.norm(self.velocity) ** 2)
        return kinetic_energy


    def calculate_momentum(self) :
        """
        Calculates and returns the momentum of a Body3D instance

        :return p: returns momentum
        """
        momentum = self.mass * self.velocity
        return momentum


    def update_position(self, dt) :
        """
        Calculates and updates the new position of a Body3D instance to 1st order

        :param dt: time-step
        """
        self.position = self.position + dt * self.velocity


    def update_2nd_position(self, dt, force) :
        """
        Calculates and updates the position of a Body3D instance to 2nd order

        :param dt: time-step
        :param force: force on body
        """
        self.position = self.position + dt * self.velocity + (dt ** 2) * (force/(2 * self.mass))


    def update_velocity(self, dt, force) :
        """
        Updates the velocity of a Body3D instance to 1st order

        :param dt: time-step
        :param force: force on body
        """
        self.velocity = self.velocity + dt * (force/self.mass)


    @staticmethod
    def new_body3d(input_file) :
        """
        Initialises a Body3D instance when given an input file handle.
        
        The input file should contain one line per body in the following format:
        label <mass> <x> <y> <z> <vx> <vy> <vz>
        
        :param input_file: Readable file handle in the above format

        :return Body3D: instance label mass position velocity
        """
        try :

            data = input_file.readline()
            lines = data.split()

            label = str(lines[0])
            mass = float(lines[1])

            x = float(lines[2])
            y = float(lines[3])
            z = float(lines[4])
            position = np.array(x, y, z)

            v_x = float(lines[5])
            v_y = float(lines[6])
            v_z = float(lines[7])
            velocity = (v_x, v_y, v_z)

            return Body3D(label, mass, position, velocity)

        except IndexError :

            print("Error: Incorrect file format.")


    @staticmethod
    def calculate_system_kinetic_energy(body3d_list) :
        """
        Returns the total kinetic energy of the system as a float

        :param body3d_list: list in which each item is a Body3D instance

        :return system_kinetic_energy: sum of each body's kinetic energy
        """

        system_kinetic_energy = 0

        for body in body3d_list :

            kinetic_energy = body.kinetic_e()
            system_kinetic_energy += kinetic_energy

        return float(system_kinetic_energy)


    @staticmethod
    def calculate_centre_of_mass_velocity(body3d_list) :
        """
        Computes the total mass and CoM velocity of a list of Body3D's

        :param body3d_list: list in which each item is a Body3D instance

        :return total_mass: The total mass of the system 
        :return com_vel: Centre-of-mass velocity
        """
        total_mass = 0
        centre_of_mass_velocity = 0
        total = 0

        for body in body3d_list :

            body_mass = body.mass
            total_mass += body_mass

        for body in body3d_list :

            body_velocity = body.vel
            mass_x_velocity = body_mass * body_velocity
            total += mass_x_velocity
            centre_of_mass_velocity = total / total_mass

        return total_mass, centre_of_mass_velocity
