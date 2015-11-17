from sys import stdin, exit
from datetime import datetime
import random
import math

print_debug = False
print_time = True


def elapsed_time_function(function_name, time0):
    time1 = datetime.now()
    elapsed_time = time1 - time0
    s = 'Elapsed time %s: %d seconds and %d milliseconds' % (function_name, elapsed_time.seconds, elapsed_time.microseconds/1000)
    if print_debug:
        print(s)
    return [elapsed_time.seconds, elapsed_time.microseconds/1000]


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
    if print_debug:
        print("Miller Rabin")
    if N == 1:
        return False
    if N & 1 == 0:
        if N == 2:
            return True
        return False
    if print_debug:
        print('N is odd')
    s = 0
    t = N-1
    while t & 1 == 0:
        t /= 2
        s += 1
    if print_debug:
        print('s = %d' % s)
        print('t = %d' % t)
    for i in range(20):
        if print_debug:
            print('Test number %d' % i)
        prob_prime = False
        a = random.randrange(1, N)
        if print_debug:
            print('Random number a = %d' % a)
        u = powermod(a, t, N)
        if print_debug:
            print('a^t mod N = %d' % u)
        if u == 1:
            prob_prime = True
        if u == N-1:
            prob_prime = True
        for j in range(s-1):
            u = u*u % N
            if print_debug:
                print('u^2 mod N = %d' % u)
            if u == N-1:
                prob_prime = True
        if print_debug:
            print('prob_prime : %d' % prob_prime)
        if not(prob_prime):
            return False
    return True


def pollard_rho(N, param):
    if print_debug:
        print("Pollard rho")
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


def brent(n):
    if n % 2 == 0:
        return 2
    x, c, m = random.randrange(0, n), random.randrange(1, n), random.randrange(1, n)
    y, r, q = x,1, 1
    g, ys = 0, 0
    while(True):
        x = y
        for i in range(r):
            y, k = (y*y+c)%n, 0
        while(True):
            ys=y
            for i in range(min(m,r-k)):
                y, q = (y*y+c)%n, q*abs(x-y) % n
            g, k = gcd(q,n), k+m
            if k>= r or g>1:
                break
        r *= 2
        if g > 1:
            break
    if g == n:
        while(True):
            ys, g = (x*x+c)%n, gcd(abs(x-ys),n)
            if g > 1:
                break
    return g, n/g


def fermat1(N):
    i = 1
    while(i <= 3):
        x = math.ceil(math.sqrt(i*N))
        j = x*x % N
        r = 0
        sqrtj = math.sqrt(j)
        while r <= 2 or math.ceil(sqrtj) == sqrtj:
            x += 1
            j = x*x % N
            sqrtj = math.sqrt(j)
            r += 1
            if math.ceil(sqrtj) == sqrtj:
                x = gcd(N, x-sqrtj)
                return (x, N/x)
        i +=1
    return (-1, -1)



def fermat2(N):
    for a in range(int(math.ceil(math.sqrt(N))), (N+9)/6):
        b = math.sqrt(a*a-N)
        if math.floor(b) == b:
            b = int(b)
            return (a-b, N/(a-b))

def fermat3(N):
    a = int(math.ceil(math.sqrt(N)))
    b = a*a-N
    sqrtb = math.sqrt(b)
    while math.ceil(sqrtb) != sqrtb:
        b += 2*a+1
        a += 1
        sqrtb = math.sqrt(b)
    return (a-int(sqrtb), a+int(sqrtb))


def factorize(N):
    if print_debug:
        print("Factorize")
    factors = []
    tofactor = []
    while N & 1 == 0:
        N /= 2
        factors.append(2)
    if print_debug:
        print('N without all the factors 2 = %d' % N)
    if N == 1:
        if print_debug:
            print('N = 1')
    elif miller_rabin(N):
        factors.append(N)
    else:
        tofactor.append(N)
        while(tofactor):
            newtofactor = []
            for f in tofactor:
                p0, p1 = brent(f)
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
         if print_debug:
            print('\ninput = %s\n' % userinput)
         v = factorize(userinput)
         s = len(v)
         if s > 0:
             current = v[0]
             exp = 0
             it = 1
             for f in v:
                 if (f==current):
                     exp += 1
                 else:
                     print ('%d^%d' % (current, exp)),
                     current = f
                     exp = 1
                 if it == s:
                     print ('%d^%d' % (current, exp))
                 it += 1
         userinput = int(stdin.readline())
    return 0

if __name__ == "__main__":
    a,b = fermat3(561)
    exit(main())
