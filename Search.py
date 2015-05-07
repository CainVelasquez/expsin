from ExpoSin import *
import math

def searchDT(expsinclass, dT, mu):
    range_ = expsinclass.tany1Range()
    a = range_[0]
    b = range_[1]
    tof_a = expsinclass.createExpoSin(a).dT(expsinclass.psi, mu)
    tof_b = expsinclass.createExpoSin(b).dT(expsinclass.psi, mu)
    c = b - (tof_b - dT) * (b - a) / (tof_b - tof_a)
    tof_c = expsinclass.createExpoSin(c).dT(expsinclass.psi, mu)

    while math.fabs(tof_c - dT) > 1.0e-6:
        if tof_a * tof_c > 0: #same sign
            a = c
        elif tof_b * tof_c > 0: #same sign
            b = c
        tof_a = expsinclass.createExpoSin(a).dT(expsinclass.psi, mu)
        tof_b = expsinclass.createExpoSin(b).dT(expsinclass.psi, mu)
        c = b - (tof_b - dT) * (b - a) / (tof_b - tof_a)
        tof_c = expsinclass.createExpoSin(c).dT(expsinclass.psi, mu)

    return expsinclass.createExpoSin(c)
