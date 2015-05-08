import math
from ExpoSin import ExpoSin

class ClassExpoSin(object):

    """
    Represents the class of sinusoids defined by S_k2[r1, r2, psi, N].

    An ExpoSin object can be constructed with this class using an
    initial tan(y1).
    """

    def __init__(self, k2, r1, r2, angle, N=0):
        self.k2 = float(k2)
        self.r1 = float(r1)
        self.r2 = float(r2)
        self.N = N
        self.psi = 2 * math.pi * N + angle

    def tany1Range(self):
        """Calculate the allowable range for tan(y1)."""
        # unpack for easy reading
        k2 = self.k2
        r1 = self.r1
        r2 = self.r2
        psi = self.psi

        logr1r2 = math.log(r1 / r2)
        cosk2O = math.cos(k2 * psi)

        delta = 2*(1-cosk2O)/k2**4 - logr1r2**2
        if delta < 0: # no feasible trajectories
            return None

        tany1min = k2/2 * (-logr1r2 / math.tan(k2*psi/2) - math.sqrt(delta))
        tany1max = k2/2 * (-logr1r2 / math.tan(k2*psi/2) + math.sqrt(delta))

        return tany1min, tany1max

    def createExpoSin(self, tany1):
        """Return a single, fully-constrained exponential sinusoid object."""
        # unpack for easy reading
        k2 = self.k2
        r1 = self.r1
        r2 = self.r2
        psi = self.psi

        logr1r2 = math.log(r1 / r2)
        sink2O = math.sin(k2 * psi)
        cosk2O = math.cos(k2 * psi)

        k1_sqr = ((logr1r2 + tany1 / k2 * sink2O)/(1 - cosk2O))**2 + (tany1 / k2)**2
        k1_sign = (logr1r2 + tany1 / k2 * sink2O)/(1 - cosk2O)
        if k1_sign < 0:
            k1 = -math.sqrt(k1_sqr)
        else:
            k1 = math.sqrt(k1_sqr)

        phi = math.acos(tany1/k1/k2)

        k0 = r1/math.exp(k1*math.sin(phi))

        return ExpoSin(k0, k1, k2, phi)
