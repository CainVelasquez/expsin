import math
from ClassExpoSin import ClassExpoSin
from ExpoSin import graph2DExpoSins, graph3DExpoSin
from PyKEP import DAY2SEC
from Vector import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D

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
        self.tany1 = solution[2]
        self.tany2 = sum(tany1Range) - self.tany1

    def impulses(self):
        """Calculate total impulse and maximum acceleration required along the arc."""
        dt = self.tof / 100
        theta = 0.0
        impulse = 0.0
        maxA = 0.0
        for i in range(100):
            requiredA = math.fabs(self.fittedExpoSin.requiredA(theta, self.mu))
            if requiredA > maxA:
                maxA = requiredA
            impulse += dt * requiredA
            theta += dt * self.fittedExpoSin.thetadot(theta, self.mu)
        return impulse, maxA

    def boundaryV(self):
        """Calculate the initial boundary velocities."""

        v1mag = self.fittedExpoSin.vmag(0.0, self.mu)
        v2mag = self.fittedExpoSin.vmag(self.psi, self.mu)

        v1_t = 1.0 / math.sqrt(1.0 + self.tany1**2) * v1mag # cosy1 * v1mag = tangential component
        v1_r = self.tany1 / math.sqrt(1.0 + self.tany1**2) * v1mag # siny1 * v1mag = radial component

        v2_t = 1.0 / math.sqrt(1.0 + self.tany2**2) * v2mag
        v2_r = self.tany2 / math.sqrt(1.0 + self.tany2**2) * v2mag

        normal = cross(self.r1, self.r2)

        dir1_t = unit(cross(normal, self.r1))
        dir1_r = unit(self.r1)

        dir2_t = unit(cross(normal, self.r2))
        dir2_r = unit(self.r2)

        v1 = add(scale(dir1_t, v1_t), scale(dir1_r, v1_r))
        v2 = add(scale(dir2_t, v2_t), scale(dir2_r, v2_r))

        if self.lw:
            v1 = scale(v1, -1.0)
            v2 = scale(v2, -1.0)

        return v1, v2

    def graph2DReduced(self, ax):
        graph2DExpoSins(ax, [self.fittedExpoSin], self.classExpoSin.psi)
        ax.plot([0.0], [self.r1mag],'ko')
        ax.plot([self.psi], [self.r2mag],'ro')
        ax.grid(True)
        ax.set_title('2D Reduced Problem')

    def graph3D(self, ax):
        graph3DExpoSin(ax, self.fittedExpoSin, self.classExpoSin.psi, self.r1, self.r2, self.lw)
        ax.plot([0], [0], [0], 'yo')
        ax.plot([self.r1[0]], [self.r1[1]], [self.r1[2]], 'ko')
        ax.plot([self.r2[0]], [self.r2[1]], [self.r2[2]], 'ro')
        #v1, v2 = self.boundaryV()
        #r1p = add(scale(v1, 5000000.0), self.r1)
        #r2p = add(scale(v2, 5000000.0), self.r2)
        #ax.plot([self.r1[0], r1p[0]], [self.r1[1], r1p[1]], [self.r1[2], r1p[2]], 'k')
        #ax.plot([self.r2[0], r2p[0]], [self.r2[1], r2p[1]], [self.r2[2], r2p[2]], 'r')
        ax.set_title('3D Trajectory')

def searchDT(expsinclass, dT, mu):
    """
    Search for a given dT using Regula Falsi method.

    Occasionally has issues with precision - the stop criterion can be relaxed,
    since the solver is usually dealing in hundreds of days while the current
    criterion looks for several-second precision! If you are getting the exception
    it throws, try this first.
    """
    range_ = expsinclass.tany1Range()
    a = range_[0]
    b = range_[1]
    tof_a = expsinclass.createExpoSin(a).dT(expsinclass.psi, mu)
    tof_b = expsinclass.createExpoSin(b).dT(expsinclass.psi, mu)

    # we have to limit the maximum TOF to ~1000 years.
    # there is a significant loss of precision in the rootfinding that causes
    # it to fail completely if it is not restricted.
    if tof_b is max(tof_b, tof_a):
        while tof_b > DAY2SEC * 365 * 1000:
            b *= 0.90
            tof_b = expsinclass.createExpoSin(b).dT(expsinclass.psi, mu)
    else:
        while tof_a > DAY2SEC * 365 * 1000:
            a *= 0.90
            tof_a = expsinclass.createExpoSin(a).dT(expsinclass.psi, mu)

    # sanity check
    if min(tof_a, tof_b) > dT or max(tof_a, tof_b) < dT:
        raise Exception('Time of flight not possible!')
    c = b - (tof_b - dT) * (b - a) / (tof_b - tof_a)
    tof_c = expsinclass.createExpoSin(c).dT(expsinclass.psi, mu)

    overflow = 10000

    while math.fabs(tof_c - dT) > 1.0 * 86400 and overflow >= 0:
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
        raise Exception('Failed to find required dT in alotted iterations. %.1f, %.1f, %.1f' % (a, b, c))

    return expsinclass.createExpoSin(c), tof_c, c
