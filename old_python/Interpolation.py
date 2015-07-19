class Polynomial(object):
    """
    Polynomial defined by c[i] x^i + ...
    """
    def __init__(self, coeff):
        self.coeff = coeff
    def __str__(self):
        string = 'y(x) = '
        for i in range(len(self.coeff) - 1, 0, -1):
            string += '({:.2e})'.format(self.coeff[i]) + '*x^%i + ' % i
        string += '{:.2e}'.format(self.coeff[0])
        return string
    def scaled(self, a):
        ret = []
        for i in range(len(self.coeff)):
            ret.append(self.coeff[i] * a)
        return Polynomial(ret)
    def multiply(self, poly):
        new_coeff_sz = len(self.coeff) + len(poly.coeff) - 1
        new_coeff = [0] * new_coeff_sz
        for i in range(len(self.coeff)):
            for j in range(len(poly.coeff)):
                power = i+j
                new_coeff[power] += self.coeff[i] * poly.coeff[j]
        return Polynomial(new_coeff)
    def add(self, poly):
        if len(self.coeff) > len(poly.coeff):
            large = self
            small = poly
        else:
            large = poly
            small = self
        new_coeff = [0] * len(large.coeff)
        for i in range( len(small.coeff) ):
            new_coeff[i] = self.coeff[i] + poly.coeff[i]
        for i in range( len(small.coeff), len(large.coeff) ):
            new_coeff[i] = large.coeff[i]
        return Polynomial(new_coeff)
    def evl(self, x):
        ret = 0.0
        x_ = 1.0
        for i in range(len(self.coeff)):
            ret += self.coeff[i] * x_
            x_ *= x
        return ret
    def integrate(self, a, b):
        ret = 0.0
        a_ = a
        for i in range(len(self.coeff)):
            ret -= self.coeff[i] * a_ / (i+1)
            a_ *= a
        b_ = b
        for i in range(len(self.coeff)):
            ret += self.coeff[i] * b_ / (i+1)
            b_ *= b
        return ret

class Interpolator(object):
    """
    Lagrangian interpolation method
    """
    def __init__(self, x_vals, y_vals):
        # lengths should be equal
        self.x_vals = x_vals
        self.y_vals = y_vals
        self.setPolynomial()
    def getY(self, x):
        return self.poly.evl(x)
    def setPolynomial(self):
        ret = Polynomial([0]) # y(x) = 0
        for i in range(len(self.x_vals)):
            ret = ret.add(self.basisPoly(i).scaled(self.y_vals[i]))
        self.poly = ret
    def basisPoly(self, j):
        p = self.numPoly(j)
        return p.scaled(1.0 / self.denomPoly(j))
    def denomPoly(self, j):
        prod = 1.0
        for i in range(len(self.x_vals)):
            if not i == j:
                prod *= self.x_vals[j] - self.x_vals[i]
        return prod
    def numPoly(self, j):
        prod = Polynomial([1]) # y(x) = 1
        for i in range(len(self.x_vals)):
            if not i == j:
                prod = prod.multiply(Polynomial([-self.x_vals[i], 1])) # y(x) = x - x_m
        return prod

def chebyshev_nodes(num, a, b):
    from math import cos, pi
    ret = []
    for i in range(1, num + 1):
        ret.append(0.5 * (a + b) + 0.5 * (b - a) * cos((2 * i - 1) * pi / 2.0 / num))
    return ret
