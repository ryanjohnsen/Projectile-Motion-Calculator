# Projectile motion solver for physics

import sys
import traceback
import math
from quadratic import quadratic_solver
from quadratic import NonrealSolutionError


# Projectile motion: a form of motion experienced by an object or
# particle (a projectile) that is projected near the Earth's surface and moves
# along a curved path under the action of gravity only (the effects of air
# resistance are assumed to be negligible). -Wikipedia

# Functionality loosely based off of
# https://www.omnicalculator.com/physics/projectile-motion

# Notes:
# horizontal acceleration = 0
# vertical acceleration = -g (gravity ~ 9.80)
# horizontal velocity is constant
# vertical velocity at max height = 0

# Variables:
# v0 -- initial velocity
# z -- launch angle (0, 90) non-inclusive for our purposes)
# y0 -- initial height (dy = -y0 if we take y1 = 0)
# t -- total time
# dx -- horizontal distance traveled (range)
# hmax -- maximum vertical height

# Equations of horizontal motion:
# vx = dx / dt
# dx = (v0 * cos(z)) * dt

# Equations of vertical motion:
# dy = (v0 * sin(z)) * dt - 0.5 * g * dt^2
# vy = v0 * sin(z) - g * dt
# vy^2 = (v0 * sin(z))^2 - 2 * g * dy
# dy = y1 - y0

# Trajectory equation:
# (y position as a function of initial speed, x position, and launch angle)
# y = tan(z) * x - [(g * x^2) / (2 * v0 * cos(z)) ^ 2]

# Horizontal range equation (if y0 = 0):
# R = (v0^2 / g) * sin(2z)

# MAIN EQUATIONS:
# (1) hmax = y0 + (v0 * sin(z))^2 / 2g
# (2) v0 * cos(z) = dx / dt
# (3) -y0 = v0 * sin(z) * t - (g/2) * t^2

class Projectile:
    def __init__(self, v0=None, angle=None, y0=None, time=None, range=None, hmax=None, gravity=9.80):
        """
            Ideally there should only be two or three arguments, but entering
            more or less is still supported. However, if more than three
            arguments (not including gravity) are given, then one or more of
            them will be overwritten.
        """
        try:
            # Initial velocity (speed)
            if v0 == None:
                self.v0 = v0
            elif v0 > 0:
                self.v0 = v0
            else:
                raise ValueError("Initial velocity must be greater than zero!")


            # Launch angle in degrees --> radians for calculations
            if angle == None:
                self.z = angle
            elif angle > 0 and angle < 90:
                self.z = math.radians(angle)
            else:
                raise ValueError("Angle is not in acceptable range (0, 90)!")

            # Initial vertical position -- determined by assuming y1 = 0
            # Change in total vertical position (dy = y1 - y0 = 0 - y0)
            if y0 == None:
                self.y0 = y0
            # Checking for TypeError
            else:
                if y0 > 0:
                    pass
                self.y0 = y0

            # Total time that projectile is in the air
            if time == None:
                self.t = time
            elif time > 0:
                self.t = time
            else:
                raise ValueError("Time must be greater than zero!")

            # Range of projectile (x1 - x0)
            if range == None:
                self.dx = range
            elif range > 0:
                self.dx = range
            else:
                raise ValueError("Range must be greater than zero!")

            # Maximum vertical height achieved by projectile
            if hmax == None:
                self.hmax = hmax
            elif hmax > 0:
                if y0 == None:
                    self.hmax = hmax
                else:
                    # Max height must be higher than y0 if launch angle > 0
                    if hmax > y0:
                        self.hmax = hmax
                    else:
                        raise ValueError("Max height can't be less than or equal to initial height!")
            else:
                # Since y1 = 0, max height must be greater than zero
                raise ValueError("Max height must be greater than zero!")

            # Gravity acting on projectile
            if gravity > 0:
                self.g = gravity
            else:
                raise ValueError("Gravity must be greater than zero!")

            # Seeing what variables we can get
            results = self.calculate_all()

            # Showing if parameters could/couldn't be determined and why
            print(results)

        except TypeError:
            # TypeError occurs when argument isn't a float or int type
            exc_type, exc_value, exc_tb = sys.exc_info()
            err_lst = traceback.format_exception(exc_type, exc_value, exc_tb)
            # Getting variable name in if statement ['if', 'var', '>', '0']
            strlst = err_lst[1].replace('\n', '').split('    ')[1].split()[1]
            err_dict = {'v0': 'Initial velocity', 'angle': 'Launch angle', 'y0': 'Initial height',
                        'time': 'Total time', 'range': 'Range', 'hmax': 'Max height', 'gravity': 'Gravity'
            }
            print(err_dict[strlst] + ' must be a numeric type!')
            sys.exit()

    def __str__(self):
        lst = [str(attribute) for attribute in self.__dict__.values()]
        info = '-----Projectile-----\nInitial Velocity (m/s): ' + lst[0] + '\n'
        info += 'Launch Angle (radians): ' + lst[1] + '\n' + 'Initial Height (m): '
        info += lst[2] + '\nTotal Time (s): ' + lst[3] + '\n' + 'Horizontal Range (m): '
        info += lst[4] + '\n' + 'Maximum Vertical Height (m): ' + lst[5]
        info += '\nGravity (m/s^2): ' + lst[6]
        return info

    def calculate_all(self):
        """Finding out the truth"""
        # Making a list of strings of the class attributes that are known
        knowns_list = [str(attribute) for attribute in self.__dict__.keys() if self.__dict__[attribute] != None]

        # Once we have three parameters, we can determine if we can solve for
        # the other parameters. If we have more than three parameters,
        # all parameters can be determined if they don't contradict each other
        if len(knowns_list) == 1 or len(knowns_list) == 2:
            # We just have gravity or one other variable
            return "No additional parameters could be calculated!"

        elif len(knowns_list) == 3:
            # We may be able to find another variable
            if 'y0' in knowns_list and 't' in knowns_list:
                # We can arrange main equations 1 and 3 to solve for max height
                # hmax = y0 + (((-y0 + (g * t^2) / 2) / (t)) ^ 2) / (2 * g)
                calculation = self.y0 + (((-self.y0 + ((self.g * self.t ** 2) / 2)) / self.t) ** 2) / (2 * self.g)
                if calculation != self.y0:
                    self.hmax = calculation
                    return "Only max height could be determined!"
                else:
                    return "Max height cannot be equal to initial height!"

            if 'y0' in knowns_list and 'hmax' in knowns_list:
                # We might be able to find time
                # A nonreal solution can only occur if hmax < y0, which is impossible
                t1, t2 = quadratic_solver(-self.g / 2, math.sqrt(2 * self.g * (self.hmax - self.y0)), self.y0)
                # We want the nonnegative solution
                if t1 > t2:
                    self.t = t1
                else:
                    self.t = t2

                return "Only total time could be determined"

            # If no values can be found, return basic message
            return "No additional parameters could be calculated"

        else:
            # We have at least three of the six important variables,
            # so we can start trying to calculate all values

            # List of all three variable trios that can calculate all parameters
            complete_vars = [('v0', 'z', 'y0'),   # [self.find_hmax(), self.find_t(2), self.find_dx()],
                             ('v0', 'z', 't'),   # [self.find_dx(), self.find_y0(2), self.find_hmax()],
                             ('v0', 'z', 'dx'),  # [self.find_t(1), self.find_y0(2), self.find_hmax()],
                             ('v0', 'z', 'hmax'), # [self.find_y0(1), self.find_t(2), self.find_dx()],
                             ('v0', 'y0', 't'),   # [self.find_z(3), self.find_dx(), self.find_hmax()],
                             ('v0', 'y0', 'hmax'), # [self.find_z(1), self.find_t(2), self.find_dx()],
                             ('v0', 't', 'dx'),   # [self.find_z(2), self.find_y0(2), self.find_hmax()],
                             ('z', 'y0', 't'),   # [self.find_v0(3), self.find_hmax(), self.find_dx()],
                             ('z', 'y0', 'hmax'), # [self.find_v0(1), self.find_t(2), self.find_dx()],
                             ('z', 't', 'dx'),   # [self.find_v0(2), self.find_y0(2), self.find_hmax()]
                             ('z', 'y0', 'dx'),  # [self.find_t(3), self.find_v0(3), self.find_hmax()]
                             ('y0', 't', 'dx'),  # [self.find_z(4), self.find_v0(2), self.find_hmax()]
                             ('y0', 'dx', 'hmax')  # [self.find_t(4), self.find_z(4), self.find_v0(1)]

                             # If we know more than three variables, then one of
                             # the trios will be a subset of those variables and
                             # thus we can solve/recalculate the three other parameters
            ]

            for i in range(len(complete_vars)):
                # Loop through each of the trios until we find suitable subset
                var_count = 0
                for var in knowns_list:
                    if var in complete_vars[i]:
                        var_count += 1
                if var_count == 3:
                    # Calling appropriate functions once we have enough
                    # parameters to solve/recalculate the other parameters
                    try:
                        # This is ugly but having functions stored in a dictionary
                        # or list is annoying since indexing or iterating over
                        # them causes all functions to be executed, so this is the
                        # best way to call the appropriate functions I could think of
                        if i == 0:
                            run_functions = [func for func in [self.find_hmax(), self.find_t(2), self.find_dx()]]
                        elif i == 1:
                            run_functions = [func for func in [self.find_dx(), self.find_y0(2), self.find_hmax()]]
                        elif i == 2:
                            run_functions = [func for func in [self.find_t(1), self.find_y0(2), self.find_hmax()]]
                        elif i == 3:
                            run_functions = [func for func in [self.find_y0(1), self.find_t(2), self.find_dx()]]
                        elif i == 4:
                            run_functions = [func for func in [self.find_z(3), self.find_dx(), self.find_hmax()]]
                        elif i == 5:
                            run_functions = [func for func in [self.find_z(1), self.find_t(2), self.find_dx()]]
                        elif i == 6:
                            run_functions = [func for func in [self.find_z(2), self.find_y0(2), self.find_hmax()]]
                        elif i == 7:
                            run_functions = [func for func in [self.find_v0(3), self.find_hmax(), self.find_dx()]]
                        elif i == 8:
                            run_functions = [func for func in [self.find_v0(1), self.find_t(2), self.find_dx()]]
                        elif i == 9:
                            run_functions = [func for func in [self.find_v0(2), self.find_y0(2), self.find_hmax()]]
                        elif i == 10:
                            run_functions = [func for func in [self.find_t(3), self.find_v0(3), self.find_hmax()]]
                        elif i == 11:
                            run_functions = [func for func in [self.find_z(4), self.find_v0(2), self.find_hmax()]]
                        else:
                            run_functions = [func for func in [self.find_t(4), self.find_z(4), self.find_v0(1)]]

                    except TypeError:
                        # If one of the methods runs into an error and the
                        # parameter can't be determined, the parameter will
                        # remain a NoneType and cause a TypeError in another method
                        return "Variables entered do not work together!"

                    # If we made it this far, then all parameters should have
                    # been calculated and thus a success statement will be printed
                    post_known_list = [str(attribute) for attribute in self.__dict__.keys() if self.__dict__[attribute] != None]
                    if len(post_known_list) == 7:
                        return "All variables have been determined!"
                    else:
                        return "Something has gone horribly wrong!"

            # If three parameters are entered and they don't match any of the trios
            return "More variables are needed to determine all parameters!"

    def find_v0(self, equation):
        if equation == 1:
            self.v0 = math.sqrt(2 * self.g * (self.hmax - self.y0)) / math.sin(self.z)

        elif equation == 2:
            self.v0 = self.dx / (self.t * math.cos(self.z))

        else:
            calculation = (-self.y0 + (self.g / 2) * self.t ** 2) / (math.sin(self.z) * self.t)
            if calculation > 0:
                self.v0 = calculation
            else:
                print("Initial velocity can't be less than or equal to zero!")

    def find_z(self, equation):
        try:
            if equation == 1:
                # arcsin domain [0, 1] --> [0, 90]
                self.z = math.asin(math.sqrt(2 * self.g * (self.hmax - self.y0)) / self.v0)

            elif equation == 2:
                # arccos domain [0, 1] --> [90, 0]
                calculation = math.acos(self.dx / (self.v0 * self.t))
                if calculation == 0:
                    raise ValueError("Launch angle cannot be zero!")
                else:
                    self.z = calculation

            elif equation == 3:
                # arcsin domain [0, 1] --> [0, 90]
                calculation = math.asin((-self.y0 + (self.g / 2) * self.t ** 2) / (self.v0 * self.t))
                if calculation == 0:
                    raise ValueError("Launch angle cannot be zero!")
                elif calculation == math.pi / 2:
                    raise ValueError("Launch angle cannot be 90 degrees!")
                else:
                    self.z = calculation

            else:
                # arctan subdomain [0, inf) --> [0, 90)
                calculation = math.atan(((self.g / 2) * self.t ** 2 - self.y0) / self.dx)
                if calculation == 0:
                    raise ValueError("Launch angle cannot be zero!")
                elif calculation < 0:
                    raise ValueError("Launch angle cannot be less than zero!")
                else:
                    self.z = calculation

        except ValueError:
            exc_type, exc_value, exc_tb = sys.exc_info()
            lst = traceback.format_exception_only(exc_type, exc_value)
            error_msg = lst[0].strip().split(': ')[1]
            if error_msg == "math domain error":
                print("Angle is not in range (0, 90)")
            else:
                print(error_msg)

    def find_y0(self, equation):
        if equation == 1:
            self.y0 = self.hmax - ((self.v0 * math.sin(self.z)) ** 2) / (2 * self.g)
        else:
            self.y0 = (self.g / 2) * self.t ** 2 - self.v0 * math.sin(self.z) * self.t

    def find_t(self, equation):
        if equation == 1:
            self.t = self.dx / (self.v0 * math.cos(self.z))

        elif equation == 2:
            if self.y0 == 0:
                # Main equation (3) simplifies to t = 2*v0*sin(z) / g
                self.t = (2 * self.v0 * math.sin(self.z)) / self.g
            else:
                try:
                    t1, t2 = quadratic_solver(-self.g / 2, self.v0 * math.sin(self.z), self.y0)
                    # We want the nonnegative/greater solution
                    if t1 > t2:
                        self.t = t1
                    else:
                        self.t = t2

                except NonrealSolutionError:
                    print('Nonreal solution')

        elif equation == 3:
            try:
                self.t = math.sqrt((self.dx * math.tan(self.z) + self.y0) / (self.g / 2))

            except ValueError:
                # dx*tan(z) - self.y0 < 0
                print("Projectile is too low")

        else:
            t1, t2 = quadratic_solver(-self.g / 2, math.sqrt(2 * self.g * (self.hmax - self.y0)), self.y0)
            # Since hmax > y0, it isn't possible to get nonreal solutions
            # We want the nonnegative/greater solution
            if t1 > t2:
                self.t = t1
            else:
                self.t = t2

    def find_dx(self):
        calculation = self.v0 * math.cos(self.z) * self.t
        if calculation > 0:
            self.dx = calculation
        else:
            print("Range cannot be less than zero!")

    def find_hmax(self):
        calculation = self.y0 + ((self.v0 * math.sin(self.z)) ** 2) / (2 * self.g)
        if calculation > 0:
            self.hmax = calculation
        else:
            print("Max height cannot be less than zero!")

# SETTERS/GETTERS
    def set_v0(self, v0):
        if v0 == None:
            self.v0 = v0
        elif type(v0) == int or type(v0) == float:
            if v0 > 0:
                self.v0 = v0
            else:
                raise ValueError("Initial velocity must be greater than zero!")
        else:
            raise TypeError("Initial velocity must be a numeric type!")


    def set_angle(self, z):
        if z == None:
            self.z = z
        elif type(z) == int or type(z) == float:
            if 0 < z < 90:
                self.z = math.radians(z)
            else:
                raise ValueError("Angle is not in acceptable range (0, 90]!")
        else:
            raise TypeError("Launch angle must be a numeric type!")

    def set_y0(self, y0):
        if y0 == None:
            self.y0 = y0
        elif type(y0) == int or type(y0) == float:
            if self.hmax == None or y0 < self.hmax:
                self.y0 = y0
            else:
                raise ValueError("Initial height cannot be greater than or equal to max height!")
        else:
            raise ValueError("Initial height must be a numeric type!")

    def set_time(self, t):
        if t == None:
            self.t = t
        elif type(t) == int or type(t) == float:
            if t > 0:
                self.t = t
            else:
                raise ValueError("Time must be greater than zero!")
        else:
            raise TypeError("Time must be a numeric type!")

    def set_range(self, dx):
        if dx == None:
            self.dx = dx
        elif type(dx) == int or type(dx) == float:
            if dx > 0:
                self.dx = dx
            else:
                raise ValueError("Range must be greater than zero!")
        else:
            raise TypeError("Range must be a numeric type!")

    def set_hmax(self, hmax):
        if hmax == None:
            self.hmax = hmax
        elif type(hmax) == int or type(hmax) == float:
            # Since y1 = 0, max height must be greater than zero
            if hmax <= 0:
                raise ValueError("Max height must be greater than zero!")
            # Max height has to be greater than initial height
            elif self.y0 == None or hmax > self.y0:
                self.hmax = hmax
            else:
                raise ValueError("Max height can't be less than or equal to initial height")
        else:
            raise TypeError("Max height must be a numeric-type or non-type!")

    def set_gravity(self, g):
        if g == None:
            self.g = g
        elif type(g) == int or type(g) == float:
            if g > 0:
                self.g = g
            else:
                raise ValueError("Gravity must be greater than zero!")
        else:
            raise TypeError("Gravity must be a numeric-type or non-type!")

    def z_degrees(self):
        if self.z != None:
            return math.degrees(self.z)

    def get_vars(self):
        """Returns a list of all the attributes"""
        return [self.v0, self.z, self.y0, self.t, self.dx, self.hmax, self.g]

def main():
    # p1 = Projectile(v0=1, angle=2, y0=-1, time=5, range=6, hmax=1, gravity=3)
    # print(p1)
    pass

if __name__ == '__main__':
    main()
