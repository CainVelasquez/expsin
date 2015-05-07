from ClassExpoSin import *
from ExpoSin import *
from ExpoSinLambert import *
from numpy import *

def testGraph():

    k2 = 0.25
    r1 = 1
    r2 = 1.5
    angle = math.pi * 2/3
    N = 1

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
    r1 = 1.0, 0.0
    r2 = 0.0, 1.5
    tof = 4.0
    mu = 1.0
    problem = ExpoSinLambert(r1, r2, tof, mu)
    graphExpoSins([problem.solexposin], problem.classExpoSin.psi)

testSolver()
