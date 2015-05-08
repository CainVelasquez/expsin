import math
from ClassExpoSin import ClassExpoSin
from ExpoSin import graphExpoSins
from PyKEP import DAY2SEC

class ExpoSinLambert(object):

    """
    Lambert solver for exponential sinusoid trajectories.

    Currently only solves the reduced 2-dimensional case.
    Currently only outputs a single solving trajectory, when in theory there
    are potentially 2N + 1 solving trajectories.
    """

    def __init__(self, r1, r2, tof, mu, N=0, lw=False, k2=0.3):
        """Currently solves entire thing within init..."""
        r1mag = math.sqrt(sum([x*x for x in r1]))
        r2mag = math.sqrt(sum([x*x for x in r2]))
        angle = math.acos(sum([x*y for x,y in zip(r1,r2)])/r1mag/r2mag)
        if lw:
            angle = math.pi * 2 - angle
        self.psi = angle + 2 * math.pi * N
        self.lw = lw
        self.tof = tof
        self.N = N
        self.k2 = k2
        self.r1 = r1
        self.r1mag = r1mag
        self.r2 = r2
        self.r2mag = r2mag
        self.mu = mu
        self.classExpoSin = ClassExpoSin(k2, r1mag, r2mag, angle, N)
        tany1Range = self.classExpoSin.tany1Range()
        if tany1Range is None:
            raise Exception('No feasible ExpoSins for given parameters')
        # We know that TOF is monotonic in tany1:
        esa = self.classExpoSin.createExpoSin(tany1Range[0]).dT(self.psi, mu)
        esb = self.classExpoSin.createExpoSin(tany1Range[1]).dT(self.psi, mu)
        min_tof = min(esa, esb)
        max_tof = max(esa, esb)
        if tof < min_tof:
            raise Exception('Time of flight too short for given parameters, min %i < tof %i < max %i (days)' % (min_tof / DAY2SEC, tof / DAY2SEC, max_tof / DAY2SEC))
        if tof > max_tof:
            raise Exception('Time of flight too long for given parameters, min %i < tof %i < max %i (days)' % (min_tof / DAY2SEC, tof / DAY2SEC, max_tof / DAY2SEC))
        solution = searchDT(self.classExpoSin, self.tof, self.mu)
        self.fittedExpoSin = solution[0]
        self.actualTOF = solution[1]

    def maxAccel(self):
        """Get thrust ceiling to follow the trajectory."""
        pass

    def totalAccel(self):
        """Calculate total impulse required along the arc."""
        pass

    def boundaryV(self):
        """Calculate the initial boundary velocities."""
        return self.fittedExpoSin.vmag(0.0, self.mu), self.fittedExpoSin.vmag(self.psi, self.mu)

    def graph(self):
        graphExpoSins([self.fittedExpoSin], self.classExpoSin.psi)

def searchDT(expsinclass, dT, mu):
    """
    Search for a given dT using Regula Falsi method.

    Occasionally has issues with precision - the stop criterion can be relaxed,
    since the solver is usually dealing in hundreds of days while the current
    criterion looks for 1-second precision! If you are getting the exception it
    throws, try this first.
    """
    range_ = expsinclass.tany1Range()
    a = range_[0]
    b = range_[1]
    tof_a = expsinclass.createExpoSin(a).dT(expsinclass.psi, mu)
    tof_b = expsinclass.createExpoSin(b).dT(expsinclass.psi, mu)
    c = b - (tof_b - dT) * (b - a) / (tof_b - tof_a)
    tof_c = expsinclass.createExpoSin(c).dT(expsinclass.psi, mu)

    overflow = 10000

    while math.fabs(tof_c - dT) > 1.0 and overflow >= 0:
        if (tof_a > 0 and tof_c > 0) or (tof_a < 0 and tof_c < 0): # same sign?
            a = c
        elif (tof_b > 0 and tof_c > 0) or (tof_b < 0 and tof_c < 0): # same sign?
            b = c
        tof_a = expsinclass.createExpoSin(a).dT(expsinclass.psi, mu)
        tof_b = expsinclass.createExpoSin(b).dT(expsinclass.psi, mu)
        c = b - (tof_b - dT) * (b - a) / (tof_b - tof_a)
        tof_c = expsinclass.createExpoSin(c).dT(expsinclass.psi, mu)
        overflow -= 1

    if overflow < 0:
        raise Exception('Failed to find required dT.')

    return expsinclass.createExpoSin(c), tof_c
