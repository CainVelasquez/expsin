'''
Plot an example transfer
'''
import PyKEP

earth = PyKEP.planet.jpl_lp('earth')
mars = PyKEP.planet.jpl_lp('mars')
k2 = 0.6
n = -1
t0 = 250.0
tof = 520.0
lw = False

r1, v1 = earth.eph(PyKEP.epoch(t0,"mjd2000"))
r2, v2 = mars.eph(PyKEP.epoch(t0+tof,"mjd2000"))

prob = PyKEP.lambert_exposin(r1, r2, tof * PyKEP.DAY2SEC, PyKEP.MU_SUN, lw, n, k2)
print prob

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

fig = plt.figure()
axis = fig.gca(projection='3d')

axis.scatter([0],[0],[0],color='y')

PyKEP.orbit_plots.plot_lambert_exposin(prob, axis)
PyKEP.orbit_plots.plot_planet(earth, PyKEP.epoch(t0,"mjd2000"), ax=axis)
PyKEP.orbit_plots.plot_planet(mars, PyKEP.epoch(t0+tof,"mjd2000"), ax=axis)

plt.show()