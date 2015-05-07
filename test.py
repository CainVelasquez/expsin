from ClassExpoSin import *
from ExpoSin import *
from ExpoSinLambert import *
from numpy import *
from PyKEP import *
from PyKEP.core import *

def testGraph():

    k2 = 0.025
    r1 = 1
    r2 = 1.5
    angle = math.pi * 3/4
    N = 5

    set = ClassExpoSin(k2, r1, r2, angle, N)
    tany1range = set.tany1Range()
    print 'tany1 range: %f, %f' % tany1range
    tanyi = linspace(tany1range[0]*0, tany1range[1], 1)
    exposins = []
    for stany in tanyi:
        exposin = set.createExpoSin(stany)
        exposins.append(exposin)
        print 'tany1: %f, dT = %f' % (stany, exposin.dT(set.psi, 1.0))

    graphExpoSins(exposins, set.psi)

def testSolver():
    kernelPath = '/home/chris/Downloads/de430.bsp'
    util.load_spice_kernel(kernelPath)
    tof = 105.0
    r1 = planet.spice('EARTH', 'SUN', 'ECLIPJ2000', 'NONE', MU_SUN, 0, 0, 0).eph(epoch(0.0))[0]
    r2 = planet.spice('VENUS', 'SUN', 'ECLIPJ2000', 'NONE', MU_SUN, 0, 0, 0).eph(epoch(tof))[0]
    problem = ExpoSinLambert(r1, r2, tof * DAY2SEC, MU_SUN, k2=0.5)
    graphExpoSins([problem.solexposin], problem.classExpoSin.psi)
    print 'dT: %f' % problem.solexposin.dT(problem.classExpoSin.psi, MU_SUN)
    print 'v1, v2: %f, %f' % problem.boundaryV()

testSolver()
