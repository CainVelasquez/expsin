'''
Check the state vector of an exposin for discrepancies
'''
import PyKEP

earth = PyKEP.planet.jpl_lp('earth')
mars = PyKEP.planet.jpl_lp('mars')
k2 = 0.6
n = 1
t0 = 250.0
tof = 520.0
lw = True

r1, v1 = earth.eph(PyKEP.epoch(t0,"mjd2000"))
r2, v2 = mars.eph(PyKEP.epoch(t0+tof,"mjd2000"))

prob = PyKEP.lambert_exposin(r1, r2, tof * PyKEP.DAY2SEC, PyKEP.MU_SUN, lw, n, k2)
print prob

'''
To look for discrepancies in propagated state, we find v by propagation as well as the implemented analytical v
'''
import random
import Vector

exps = prob.get_exposins()[0]

# any random progress into the trajectory
rand_psi = random.uniform(0.0, exps.get_psi())
r, v, a = exps.get_state(rand_psi, PyKEP.MU_SUN)

# we take the angular variation of r using given v
dpsi = Vector.mag(Vector.scale(Vector.cross(r, v), 1.0/Vector.mag(r)**2))

# and pull the next state
rp, vp, ap = exps.get_state(rand_psi+dpsi, PyKEP.MU_SUN)

# we thus have r -> rp over dt = 1.0, so we can directly extract an approximate v
v_approx = Vector.sub(rp, r)

print 'Similarity, v -> v_approx: %0.8f' % (Vector.dot(v_approx, v) / Vector.mag(v_approx) / Vector.mag(v))
print 'Similarity, vp -> v_approx: %0.8f' % (Vector.dot(v_approx, vp) / Vector.mag(v_approx) / Vector.mag(vp))
