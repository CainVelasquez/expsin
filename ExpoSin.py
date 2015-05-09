import math
import matplotlib.pyplot as plt
from Vector import *

class ExpoSin(object):

    """Represents a single, fully constrained exponential sinusoid."""

    def __init__(self, k0, k1, k2, phi):
        self.k0 = k0
        self.k1 = k1
        self.k2 = k2
        self.phi = phi

    def __str__(self):
        """Pretty print the parameters for the given exponential sinusoid."""
        return 'k0 = %f, k1 = %f, k2 = %f, phi = %f' % (self.k0, self.k1, self.k2, self.phi)

    def tany(self, theta):
        """Get the flight path angle stored as its tangent."""
        return self.k1 * self.k2 * math.cos(self.k2 * theta + self.phi)

    def r(self, theta):
        """Get the radius."""
        return self.k0 * math.exp(self.k1 * math.sin(self.k2 * theta + self.phi))

    def rdot(self, theta, mu):
        """Get the radius derivative w.r.t. time."""
        return self.r(theta) * self.k1 * math.cos(self.k2 * theta + self.phi) * self.k2 * self.thetadot(theta, mu)

    def thetadot(self, theta, mu):
        """Get the angular rate w.r.t. time."""
        return math.sqrt(mu / self.r(theta)**3 / (self.tany(theta)**2
        + self.k1 * self.k2**2 * math.sin(self.k2 * theta + self.phi) + 1))

    def dT(self, psi, mu):
        """Calculate the time to traverse the given angle."""
        #Performs a placeholder midpoint algorithm to evaluate the quadrature.
        dtheta = psi / 100
        thetai = [dtheta * t + dtheta / 2 for t in range(100)]
        integrandi = [1.0/self.thetadot(t, mu) * dtheta for t in thetai]
        integral = sum(integrandi)
        return integral

    def vmag(self, theta, mu):
        """Calculate velocity magnitude."""
        t = self.thetadot(theta, mu) * self.r(theta)
        r = self.rdot(theta, mu)
        v = math.sqrt(t**2 + r**2)
        return v

    def localA(self, theta):
        """Calculate the local acceleration required as a fraction of local g."""
        tany = self.tany(theta)
        s = math.sin(self.k2 * theta + self.phi)
        cosy = 1 / math.sqrt(1 + tany**2)
        tan2yk1k22s1 = tany**2 + self.k1 * self.k2**2 * s + 1
        return tany / 2.0 / cosy * (1 / tan2yk1k22s1 - (self.k2**2 * (1 - 2 * self.k1 * s))/tan2yk1k22s1**2)

    def requiredA(self, theta, mu):
        """Calculate the required acceleration to follow the curve."""
        return mu / self.r(theta)**2 * self.localA(theta)

def graph2DExpoSins(ax, exposins, psi):
    """Graph an array of exponential sinusoid objects."""
    thetai = [psi / 400 * t for t in range(401)]
    for exposin in exposins:
        radiusi = [exposin.r(t) for t in thetai]
        ax.plot(thetai, radiusi, 'k')

def graph3DExpoSin(ax, exposin, psi, r1, r2, backward):
    normal = cross(r1, r2)
    thetai = [psi / 400 * t for t in range(401)]
    radiusi = [exposin.r(t) for t in thetai]
    if backward:
        rti = [-r * math.sin(t) for r, t in zip(radiusi, thetai)]
    else:
        rti = [r * math.sin(t) for r, t in zip(radiusi, thetai)]
    rri = [r * math.cos(t) for r, t in zip(radiusi, thetai)]
    coordsi = [transform(rt, rr, r1, normal) for rt, rr in zip(rti, rri)]
    x, y, z = zip(*coordsi)
    ax.plot(x, y, z, 'k')
