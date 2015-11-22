from sys import stdin, exit
from math import ceil, floor, sqrt
from random import randrange


def gcd(a, b):
    while b != 0:
        tmp = b
        b = a % b
        a = tmp
    return a


def powermod(x, n, m):
    if n == 0:
        return 1
    if n == 1:
        return x
    y = powermod(x, n/2, m)
    if n & 1 == 0:
        return (y*y)%m
    else:
        return (((y*y)%m)*x)%m


def miller_rabin(N):
    if N == 1:
        return False
    if N & 1 == 0:
        if N == 2:
            return True
        return False
    s = 0
    t = N-1
    while t & 1 == 0:
        t /= 2
        s += 1
    for i in range(20):
        prob_prime = False
        a = randrange(1, N)
        u = powermod(a, t, N)
        if u == 1:
            prob_prime = True
        if u == N-1:
            prob_prime = True
        for j in range(s-1):
            u = u*u % N
            if u == N-1:
                prob_prime = True
        if not(prob_prime):
            return False
    return True


def pollard_rho(N, param):
    x = 2
    saved_x = 2
    y = 2
    saved_y = 2
    r = 1
    d = 1
    e = 0
    while d == 1:
        x = (x*x + param) % N
        temp = (y*y + param) % N
        y = (temp*temp + param) % N
        if x > y:
            r = (r*(x-y)) % N
        else:
            r = (r*(y-x)) % N
        if e % 100 == 0:
            d = gcd(r, N)
            if d == 1:
                saved_x = x
                saved_y = y
        e += 1
    if d == N:
        x = saved_x
        y = saved_y
        d = 1
        while d == 1:
            x = (x*x + param) % N
            temp = (y*y + param) % N
            y = (temp*temp + param) % N
            if x > y:
                d = gcd(x-y, N)
            else:
                d = gcd(y-x, N)
        if d == N:
            return pollard_rho(N, param + 1)
    return d, N/d


def brent(n, c):
    if n % 2 == 0:
        return 2, n/2
    m = 100
    x = 2
    y = x
    r = 5
    q = 1
    g = 0
    ys = 0
    while True:
        x = y
        k = 0
        for i in range(r):
            y = (y*y+c) % n
        while True:
            ys = y
            for i in range(min(m, r-k)):
                q = q*abs(x-y) % n
                y = (y*y+c) % n
            g = gcd(q, n)
            k += m
            if k >= r or g > 1:
                break
        r *= 2
        if g > 1:
            break
    if g == n:
        while True:
            g = gcd(abs(x-ys),n)
            ys = (ys*ys+c) % n
            if g > 1:
                break
        if g == n:
            return brent(n, c+1)
    return g, n/g


def factorize(N):
    factors = []
    tofactor = []
    while N & 1 == 0:
        N /= 2
        factors.append(2)
    if N != 1:
        if miller_rabin(N):
            factors.append(N)
        else:
            tofactor.append(N)
            while(tofactor):
                newtofactor = []
                for f in tofactor:
                    p0, p1 = brent(f, 1)
                    if miller_rabin(p0):
                        factors.append(p0)
                    else:
                        newtofactor.append(p0)
                    if miller_rabin(p1):
                        factors.append(p1)
                    else:
                        newtofactor.append(p1)
                tofactor = newtofactor
    factors.sort()
    return factors


def main():
    userinput = int(stdin.readline())
    while userinput!=0:
         v = factorize(userinput)
         s = len(v)
         if s > 0:
             current = v[0]
             exp = 0
             it = 1
             for f in v:
                 if f == current:
                     exp += 1
                 else:
                     print('%d^%d' % (current, exp)),
                     current = f
                     exp = 1
                 if it == s:
                     print('%d^%d' % (current, exp))
                 it += 1
         userinput = int(stdin.readline())
    return 0


if __name__ == "__main__":
    exit(main())
