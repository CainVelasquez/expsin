'''
This code will produce a porkchop plot of the J function over a range of launch times and flight times
'''
import math
import numpy
import PyKEP as pk
import matplotlib.pyplot as plt
import Vector

p1 = pk.planet.jpl_lp('earth')
p2 = pk.planet.mpcorb('99942   19.2   0.15 K107N 202.49545  126.41859  204.43202    3.33173  0.1911104  1.11267324   0.9223398  1 MPO164109  1397   2 2004-2008 0.40 M-v 3Eh MPCAPO     C802  (99942) Apophis            20080109')
k2 = 0.6
n = 0
isp_chem = 350
isp_lt = 3000

t0_range = [0, 5000]
tof_range = [100, 900]

@numpy.vectorize
def J(t0, tof):
    ep = pk.epoch(t0 + tof)
    r1, v1 = p1.eph(ep)
    r2, v2 = p2.eph(ep)
    resJ = []
    for lw in [False, True]:
        prob = pk.lambert_exposin(r1, r2, tof * pk.DAY2SEC, pk.MU_SUN, lw, n, k2)
        for i in range(prob.num_solutions()):
            exps = prob.get_exposins()[i]
            dv1 = Vector.mag(Vector.sub(prob.get_v1()[i], v1))
            dv2 = Vector.mag(Vector.sub(prob.get_v2()[i], v2))
            dvlt = exps.get_delta_v(pk.MU_SUN)
            resJ.append(1.0 - math.exp(-(dv1 + dv2) / 9.81 / isp_chem - dvlt / 9.81 / isp_lt))
    if len(resJ) == 0:
        return numpy.nan
    else:
        return numpy.nanmin(resJ)

t0 = numpy.linspace(t0_range[0], t0_range[1], 200)
tof = numpy.linspace(tof_range[0], tof_range[1], 200)
X, Y = numpy.meshgrid(t0, tof)
Z = J(X, Y)
mZ = numpy.ma.array(Z, mask=numpy.isnan(Z))
plt.figure()
plt.pcolor(X, Y, mZ, vmin=numpy.nanmin(Z), vmax=numpy.nanmax(Z))
plt.colorbar()
plt.show()
