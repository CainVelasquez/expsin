from ExpoSinLambert import *
r1 = [149.0e9, 0.0, 0.0]
r2 = [200.0e9, -100.0e9, 0.0]
tof = 86400 * 565.0;
mu = 1.32e20;
lw = True;
N = 0;
k2 = 0.175;
x = ExpoSinLambert(r1, r2, tof, mu, N, lw, k2)
print x.tany1
