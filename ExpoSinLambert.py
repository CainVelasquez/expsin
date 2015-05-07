import math
from ClassExpoSin import *
from ExpoSin import *
from Search import *

class ExpoSinLambert(object):
    """
    Provides a full description of a trajectory based on an exponential sinusoid
    using the parameterization of a Lambert problem.
    """
    def __init__(self, r1, r2, tof, mu, N=0, lw=False, k2=0.125):
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
        min_tof_exposin = self.classExpoSin.createExpoSin(tany1Range[0])
        max_tof_exposin = self.classExpoSin.createExpoSin(tany1Range[1])
        if tof < min_tof_exposin.dT(self.psi, mu) or tof > max_tof_exposin.dT(self.psi, mu):
            raise Exception('Time of flight not achievable with given parameters')
        self.solexposin = searchDT(self.classExpoSin, self.tof, self.mu)

    def maxAccel(self):
        pass

    def boundaryV(self):
        
        pass
