"""
Test the Lambert solver for the simplified problem. Currently capable of extracting boundary velocities.

This is paired with the proof-of-concept IPython notebook.

THIS IS NOT UP TO DATE!... See the IPython notebook instead!
"""

from ExpoSinLambert import ExpoSinLambert
from PyKEP import MU_SUN, DAY2SEC, AU
from PyKEP.core import epoch
from PyKEP.planet import jpl_lp
import math

# Simple earth to mars transfer problem
start = 0.0 # days
tof = 700.0 # days
p1 = jpl_lp('earth')
p2 = jpl_lp('mars')
N = 1
lw = False

# k2 - the free variable that describes how quickly the radius cycles with time.
# the solver in its current form is rather finnicky about your choice of k2;
# you will want to explore potentially many values (many will result in errors).
k2 = 0.5

# Solution of the ExpoSin Lambert problem
r1, v1 = p1.eph(epoch(start))
r2, v2 = p2.eph(epoch(start + tof))

problem = ExpoSinLambert(r1, r2, tof*DAY2SEC, MU_SUN, N=N, k2=k2, lw=lw)

print 'radius of origin (day 0) = %.5f AU' % (math.sqrt(sum([x*x for x in r1])) / AU)
print 'radius of target (day %i) = %.5f AU' % (tof+start, math.sqrt(sum([x*x for x in r2])) / AU)
print 'traversed angle = %i degrees, %.2f radians' % (problem.classExpoSin.psi * 180 / math.pi, problem.classExpoSin.psi)
print 'calculated dT = %.2f days' % (problem.actualTOF / DAY2SEC)
print 'solving ExpoSin: ' + str(problem.fittedExpoSin)
esv1, esv2 = problem.boundaryV()
print 'origin v1, target v2 = %.3f, %.3f km/s' % (math.sqrt(sum([x*x for x in v1])) / 1000.0, math.sqrt(sum([x*x for x in v2])) / 1000.0)
print 'esv1, esv2 = %.3f, %.3f km/s' % (esv1 / 1000.0, esv2 / 1000.0)
print 'total impulse = %.2f N-s, max accel = %.5f m/s^2' % problem.impulses()

problem.graph()
