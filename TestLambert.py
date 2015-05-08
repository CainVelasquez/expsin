"""
Test the Lambert solver for the simplified problem. Currently capable to extracting boundary velocities.
"""

from ExpoSinLambert import ExpoSinLambert
from PyKEP import MU_SUN, DAY2SEC, AU
from PyKEP.core import epoch
from PyKEP.planet import jpl_lp
import math

# Simple earth to mars transfer problem
start = 100.0 # days
tof = 900.0 # days
p1 = jpl_lp('earth')
p2 = jpl_lp('mars')
N = 1
lw = False

# k2 - the free variable that describes how quickly the radius cycles with time.
# the solver in its current form is rather finnicky about your choice of k2;
# you will want to explore potentially many values (many will result in errors).
# It can be shown that all possible k2 values are within the range
#       k2^4 = (4 / ln^2(r1/r2))
# from (11) in Izzo's paper.
k2 = 0.75

# Solution of the ExpoSin Lambert problem
r1, v1 = p1.eph(epoch(start))
r2, v2 = p2.eph(epoch(start + tof))

problem = ExpoSinLambert(r1, r2, tof*DAY2SEC, MU_SUN, N=N, k2=k2, lw=lw)

print 'radius of origin (day 0) = %f AU' % (math.sqrt(sum([x*x for x in r1])) / AU)
print 'radius of target (day %i) = %f AU' % (tof+start, math.sqrt(sum([x*x for x in r2])) / AU)
print 'traversed angle = %f degrees' % (problem.classExpoSin.psi * 180 / math.pi)
print 'solving ExpoSin: ' + str(problem.fittedExpoSin)
esv1, esv2 = problem.boundaryV()
print 'origin v1, target v2 = %f, %f km/s' % (math.sqrt(sum([x*x for x in v1])) / 1000.0, math.sqrt(sum([x*x for x in v2])) / 1000.0)
print 'esv1, esv2 = %f, %f km/s' % (esv1 / 1000.0, esv2 / 1000.0)

problem.graph()
