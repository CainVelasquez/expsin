import math

class QuadraticRoots(object):
    """
    Finds roots of a quadratic equation.
    """
    @classmethod
    def getQE(self, a, b, c):
        """
        Simple application of quadratic equation.
        """
        return (-b+math.sqrt(b**2-4*a*c))/2/a, (-b-math.sqrt(b**2-4*a*c))/2/a
    @classmethod
    def getRQE(self, a, b, c):
        """
        A recast, stabler alternative to the quadratic equation
        for noisy inputs.
        """
        a = float(a)
        b = float(b)
        c = float(c)
        if b <= 0:
            x1 = (-b + math.sqrt(b**2 - 4*a*c))/2/a
            x2 = c / a / x1
        else:
            x2 = (-b - math.sqrt(b**2 - 4*a*c))/2/a
            x1 = c / a / x2
        return x1, x2
    @classmethod
    def getJT(self, a, b, c):
        """
        Jenkins-Traub rootfinding algorithm.
        """
        xi1 = 0.0
        while math.fabs(a*xi1**2 + b*xi1 + c) > 1.0e-9:
            y = a*xi1**2 + b*xi1 + c
            dy = 2*a*xi1 + b
            xi1 += -y/dy
        xi2 = 0.0
        while math.fabs((a*xi2**2 + b*xi2 + c) / (xi2 - xi1)) > 1.0e-9:
            y = (a*xi2**2 + b*xi2 + c) / (xi2 - xi1)
            dy = ((xi2-xi1) * (2*a*xi2+b) - (a*xi2**2 + b*xi2 + c))/(xi2-xi1)**2
            xi2 += -y/dy
        return xi1, xi2
