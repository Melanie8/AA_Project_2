# lenstra's algorithm

from random import randint
from fractions import gcd

def primes(n):
    b, p, ps = [True] * (n+1), 2, []
    for p in xrange(2, n+1):
        if b[p]:
            ps.append(p)
            for i in xrange(p, n+1, p):
                b[i] = False
    return ps

def bezout(a, b):
    if b == 0: return 1, 0, a
    q, r = divmod(a, b)
    x, y, g = bezout(b, r)
    return y, x-q*y, g

def add(p, q, a, b, m):
    if p[2] == 0: return q
    if q[2] == 0: return p
    if p[0] == q[0]:
        if (p[1] + q[1]) % m == 0:
            return 0, 1, 0 # infinity
        n = (3 * p[0] * p[0] + a) % m
        d = (2 * p[1]) % m
    else:
        n = (q[1] - p[1]) % m
        d = (q[0] - p[0]) % m
    x, y, g = bezout(d, m)
    if g > 1: return 0, 0, d # failure
    z = (n*x*n*x - p[0] - q[0]) % m
    return z, (n * x * (p[0] - z) - p[1]) % m, 1

def mul(k, p, a, b, m):
    r = (0,1,0)
    while k > 0:
        if p[2] > 1:
            return p
        if k % 2 == 1:
            r = add(p, r, a, b, m)
        k = k // 2
        p = add(p, p, a, b, m)
    return r

def lenstra(n, limit):
    g = n
    while g == n:
        q = randint(0, n-1), randint(0, n-1), 1
        a = randint(0, n-1)
        b = (q[1]*q[1] - q[0]*q[0]*q[0] - a*q[0]) % n
        g = gcd(4*a*a*a + 27*b*b, n)
    if g > 1: return g # lucky factor
    for p in primes(limit):
        pp = p
        while pp < limit:
            q = mul(p, q, a, b, n)
            if q[2] > 1:
                return gcd(q[2], n)
            pp = p * pp
    return False

if __name__ == "__main__":
    print lenstra(36, 1000)
