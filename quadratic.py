# Quadratic equation solver function
# Only returns a tuple for real solutions -- nonreal solutions raise an error

import math

class NonrealSolutionError(Exception):
    """Occurs when a solution is nonreal"""
    # def __str__(self):
    #     return 'Quadratic has nonreal solutions'

def quadratic_solver(a, b, c):
    try:
        """Returns the solutions as a tuple"""
        # Quadratic equation: ax^2 + bx + c = 0
        # Quadratic formula: x = (-b +- sqrt(b^2 - 4ac)) / 2a --> where a =/= 0
        # Discriminant: b^2 - 4ac
        discriminant = b ** 2 - 4 * a * c

        # If discriminant < 0, the solutions have an imaginary component
        if discriminant < 0:
            raise NonrealSolutionError

        elif discriminant == 0:
            # The roots are equivalent
            root1 = -b / (2 * a)
            root2 = root1
        else:
            # General quadratic formula used
            root1 = (-b + math.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
            root2 = (-b - math.sqrt(b ** 2 - 4 * a * c)) / (2 * a)

        return root1, root2

    except TypeError:
        print('Arguments must be numeric types')

def main():
    #print(quadratic_solver(1/16, 9.8, -245))
    #print(quadratic_solver(5, -10, -4))
    pass

if __name__ == '__main__':
    main()
