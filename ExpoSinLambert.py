import math
from ClassExpoSin import ClassExpoSin
from ExpoSin import graph2DExpoSin, graph3DExpoSin
from PyKEP import DAY2SEC
from Vector import *

class ExpoSinLambert(object):

    """
    Lambert solver for exponential sinusoid trajectories.

    Currently only outputs a single solving trajectory, when in theory there
    are potentially 2N + 1 solving trajectories.
    """

    def __init__(self, r1, r2, tof, mu, N=0, lw=False, k2=0.3):
        """
        Set up the problem with the given Lambert problem parameterization.
        r1 := starting position (array3)
        r2 := ending position (array3)
        tof := time in seconds for the transfer
        mu := gravitational parameter
        N := number of revolutions to force before intercept
        lw := traverse long way instead of shortest angle
        k2 := free winding parameter of exponential sinusoid, choose carefully!
        """
        # Currently solves the entire thing in init...
        r1mag = mag(r1)
        r2mag = mag(r2)
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
            raise Exception('No feasible exponential sinusoid trajectories for given problem!')
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
        ABSCISSAS = 100
        dt = self.tof / ABSCISSAS
        theta = 0.0
        impulse = 0.0
        maxA = 0.0
        for i in range(ABSCISSAS):
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

        v1 = transform(v1_t, v1_r, self.r1, normal)
        v2 = transform(v2_t, v2_r, self.r2, normal)

        if self.lw:
            v1 = scale(v1, -1.0)
            v2 = scale(v2, -1.0)

        return v1, v2

    def graph2DReduced(self, ax):
        """Graph the Lambert arc for the 2D reduced case."""
        graph2DExpoSin(ax, self.fittedExpoSin, self.classExpoSin.psi)
        ax.plot([0.0], [self.r1mag],'ko')
        ax.plot([self.psi], [self.r2mag],'ro')
        ax.grid(True)
        ax.set_title('2D Reduced Problem')

    def graph3D(self, ax):
        """Graph the Lambert arc in the full 3D representation."""
        graph3DExpoSin(ax, self.fittedExpoSin, self.classExpoSin.psi, self.r1, self.r2, self.lw)
        ax.set_title('3D Trajectory')

def searchDT(expsinclass, dT, mu):
    """
    Search for a given dT using Regula Falsi method.

    Occasionally has issues with precision - the stop criterion can be relaxed,
    since the solver is usually dealing in hundreds of days while the current
    criterion looks for several-second precision! If you are getting the exception
    it throws, try this first.

    Returns the solving ExpoSin, its TOF, and its tany1
    """
    STOP_CRITERION = 1.0 # seconds
    MAX_ITERS = 10000

    range_ = expsinclass.tany1Range()
    a = range_[0]
    b = range_[1]
    tof_a = expsinclass.createExpoSin(a).dT(expsinclass.psi, mu)
    tof_b = expsinclass.createExpoSin(b).dT(expsinclass.psi, mu)

    # we have to limit the maximum TOF to ~1000 years.
    # there is a significant loss of precision in the= rootfinding that causes
    # it to fail completely if it is not restricted.
    MAX_DT = 1000 * 365.0 * DAY2SEC
    if dT > MAX_DT:
        raise Exception('Cannot search for dT over %i years, try changing the limit' % (MAX_DT / 365 / DAY2SEC))
    # Clamp the tany1 range to fit the restriction
    while tof_b > MAX_DT:
        b *= 0.95
        tof_b = expsinclass.createExpoSin(b).dT(expsinclass.psi, mu)
    while tof_a > MAX_DT:
        a *= 0.95
        tof_a = expsinclass.createExpoSin(a).dT(expsinclass.psi, mu)

    c = b - (tof_b - dT) * (b - a) / (tof_b - tof_a)
    tof_c = expsinclass.createExpoSin(c).dT(expsinclass.psi, mu)

    while math.fabs(tof_c - dT) > STOP_CRITERION and MAX_ITERS >= 0:
        if (tof_a > 0 and tof_c > 0) or (tof_a < 0 and tof_c < 0): # same sign?
            a = c
        elif (tof_b > 0 and tof_c > 0) or (tof_b < 0 and tof_c < 0): # same sign?
            b = c
        tof_a = expsinclass.createExpoSin(a).dT(expsinclass.psi, mu)
        tof_b = expsinclass.createExpoSin(b).dT(expsinclass.psi, mu)
        c = b - (tof_b - dT) * (b - a) / (tof_b - tof_a)
        tof_c = expsinclass.createExpoSin(c).dT(expsinclass.psi, mu)
        MAX_ITERS -= 1

    if MAX_ITERS < 0:
        raise Exception('Failed to find required dT in allotted iterations.')

    return expsinclass.createExpoSin(c), tof_c, c
