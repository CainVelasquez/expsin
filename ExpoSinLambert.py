import math
from ClassExpoSin import *
from ExpoSin import *

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
        min_tof = min_tof_exposin.dT(self.psi, mu)
        max_tof = max_tof_exposin.dT(self.psi, mu)
        if tof < min_tof or tof > max_tof:
            raise Exception('Time of flight not achievable with given parameters, %f < tof < %f' % (min_tof,max_tof))
        self.solexposin = searchDT(self.classExpoSin, self.tof, self.mu)

    def maxAccel(self):
        pass

    def boundaryV(self):
        # at start:
        t = self.solexposin.thetadot(0.0, self.mu) * self.solexposin.r(0.0)
        r = self.solexposin.rdot(0.0)
        v1 = math.sqrt(t**2 + r**2)
        # at end:
        t = self.solexposin.thetadot(self.psi, self.mu) * self.solexposin.r(self.psi)
        r = self.solexposin.rdot(self.psi)
        v2 = math.sqrt(t**2 + r**2)

        return v1, v2

def searchDT(expsinclass, dT, mu):
    range_ = expsinclass.tany1Range()
    a = range_[0]
    b = range_[1]
    tof_a = expsinclass.createExpoSin(a).dT(expsinclass.psi, mu)
    tof_b = expsinclass.createExpoSin(b).dT(expsinclass.psi, mu)
    c = b - (tof_b - dT) * (b - a) / (tof_b - tof_a)
    tof_c = expsinclass.createExpoSin(c).dT(expsinclass.psi, mu)

    while math.fabs(tof_c - dT) > 1.0e-6:
        if tof_a * tof_c > 0: #same sign
            a = c
        elif tof_b * tof_c > 0: #same sign
            b = c
        tof_a = expsinclass.createExpoSin(a).dT(expsinclass.psi, mu)
        tof_b = expsinclass.createExpoSin(b).dT(expsinclass.psi, mu)
        c = b - (tof_b - dT) * (b - a) / (tof_b - tof_a)
        tof_c = expsinclass.createExpoSin(c).dT(expsinclass.psi, mu)

    return expsinclass.createExpoSin(c)
