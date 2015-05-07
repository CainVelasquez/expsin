import math

class ExpoSin(object):
    """
    Represents a single exponential sinusoid.
    """
    def __init__(self, k0, k1, k2, phi):
        self.k0 = k0
        self.k1 = k1
        self.k2 = k2
        self.phi = phi

    def __str__(self):
        return 'k0= %f, k1= %f, k2= %f, phi= %f' % (self.k0, self.k1, self.k2, self.phi)

    def tany(self, theta):
        return self.k1 * self.k2 * math.cos(self.k2 * theta + self.phi)

    def r(self, theta):
        return self.k0 * math.exp(self.k1 * math.sin(self.k2 * theta + self.phi))

    def thetadot(self, theta, mu):
        return math.sqrt(mu / self.r(theta)**3 / (self.tany(theta)**2
        + self.k1 * self.k2**2 * math.sin(self.k2 * theta + self.phi) + 1))

    def dT(self, psi, mu):
        """
        Performs a placeholder midpoint algorithm to evaluate the quadrature.
        """
        dtheta = psi / 100
        thetai = [dtheta * t + dtheta / 2 for t in range(100)]
        integrandi = [1/self.thetadot(t, mu) * dtheta for t in thetai]
        integral = sum(integrandi)
        return integral

def graphExpoSins(exposins, psi):
    import matplotlib.pyplot as plt
    thetai = [psi / 400 * t for t in range(401)]
    for exposin in exposins:
        radiusi = [exposin.r(t) for t in thetai]
        ax = plt.subplot(111, polar=True)
        ax.plot(thetai, radiusi, 'k')
    ax.grid(True)
    plt.show()
