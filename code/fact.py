from sys import stdin, exit
from datetime import datetime
import random
import math
import numpy as np
from sets import Set

print_debug = True
print_time = False
first_factorisations_mode = True
max_dec = 1
first_fact = ()

def elapsed_time_function(function_name, time0):
    time1 = datetime.now()
    elapsed_time = time1 - time0
    s = 'Elapsed time %s: %d seconds and %d milliseconds' % (function_name, elapsed_time.seconds, elapsed_time.microseconds/1000)
    if print_time:
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


def brent(n, c):
    if n % 2 == 0:
        return 2, n/2
    #x, c, m = random.randrange(0, n), random.randrange(1, n), random.randrange(1, n)
    m = 100
    x = 2
    y = x
    r = 1
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


def fermat1(N):
    i = 1
    while i <= 3:
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
        i += 1
    return -1, -1


def fermat2(N):
    for a in range(int(math.ceil(math.sqrt(N))), (N+9)/6):
        b = math.sqrt(a*a-N)
        if math.floor(b) == b:
            b = int(b)
            return a-b, N/(a-b)


def fermat3(N):
    a = int(math.ceil(math.sqrt(N)))
    b = a*a-N
    sqrtb = math.sqrt(b)
    while math.ceil(sqrtb) != sqrtb:
        b += 2*a+1
        a += 1
        sqrtb = math.sqrt(b)
    return a-int(sqrtb), a+int(sqrtb)


def decompose(b, baseSize, base):
    dec = np.zeros(baseSize, dtype=np.int)
    dec_parity = np.zeros(baseSize, dtype=np.int)
    i = 0
    for f in base:
        while b % f == 0:
            dec[i] += 1
            b /= f
        dec_parity[i] = dec[i] & 1
        i += 1
    if b > 1:
        return None, None
    else:
        return dec, dec_parity


def gauss_jordan(A, baseSize):
    comb = [[i] for i in range(baseSize+1)]
    index = [i for i in range(baseSize+1)]
    r = -1
    for j in range(baseSize):
        if print_debug:
            print("A =")
            print A
            print("comb =")
            print comb
        # Find pivot
        for k in range(r+1, baseSize+1):
            if A[k, j] == 1:
                break
        if print_debug:
            print("The pivot is in position (%d, %d)" % (k, j))
        # Subtract
        if A[k, j] != 0:
             r += 1
             if k != r:
                 if print_debug:
                     print("k = %d, r = %d" % (k, r))
                 temp = A[r, :]
                 if print_debug:
                     print A[r, :]
                 A[r, :] = A[k, :]
                 if print_debug:
                     print A[r, :]
                 A[k, :] = temp
                 if print_debug:
                     print A[k, :]
                 temp = index[k]
                 index[k] = index[r]
                 index[r] = temp
                 if print_debug:
                    print("A =")
                    print A
             for i in range(baseSize+1):
                 if i != r:
                     if A[i, j] != 0:
                         for c in comb[index[r]]:
                            comb[index[i]].append(c)
                         for l in range(baseSize):
                             A[i, l] = (A[i, l] - A[r, l]) & 1
                             if print_debug:
                                 print("entry (%d, %d) becomes %d" %(i, l, A[i, l] - A[r, j]))
                         '''if not np.any(A[i, :] % 2):
                             comb[index[i]].append(index[i])
                             return comb[index[i]]'''
    if print_debug:
        print("A =")
        print A
        print("comb =")
        print comb
        print ("index =")
        print index
    return comb[baseSize]

def quadratic_sieve(N):
     base = [2, 7, 13]
     baseSize = 3
     dec = np.zeros((baseSize+1, baseSize), dtype=np.int)
     dec_parity = np.zeros((baseSize+1, baseSize), dtype=np.int)
     a = int(math.ceil(math.sqrt(N)))
     b = a*a-N
     counter = 0
     while counter <= baseSize:
         print a - 42
         d, d_parity = decompose(b, baseSize, base)
         if d is not None:
             print d
             if not np.any(d_parity % 2):
                 g = gcd(a+int(math.sqrt(b)), N)
                 return g, N/g
             dec[counter] = d
             dec_parity[counter] = d_parity
             counter += 1
         a += 1
         b = a*a-N
     # Find even combination
     comb = gauss_jordan(dec_parity, baseSize)
     d = np.zeros(baseSize, dtype=np.int)
     for x in comb:
         d[:] += dec[x,:]
     a = 1
     i = 0
     for b in base:
         a *= math.pow(b, d[i]/2)
         i += 1
     g = gcd(a-N, N)
     return g, N/g

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
                if not(first_factorisations_mode) and f <= max_dec:
                    for g in first_fact[f]:
                        factors.append(g)
                else:
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
    if first_factorisations_mode:
        return factors
    factors.sort()
    return factors


def main():
    time0 = datetime.now()
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
    elapsed_time_function("main", time0)
    return 0


def first_factorisations():
    file = open("first_factorisations.txt", "w")
    file.write("first_fact = ([], [], [2], ")
    for i in range(3, max_dec+1):
        file.write("[")
        v = factorize(i)
        for j in range(len(v)):
            if j == 0:
                file.write("%d" % v[j])
            else:
                file.write(", %d" % v[j])
        if i == max_dec:
            file.write("])")
        else:
            file.write("], ")
    file.close()

if __name__ == "__main__":
    #x = 834
    # 836 -> 0.26
    # 834 -> 0.25
    # 835 -> 0.32
    #random.seed(x)
    #if first_factorisations_mode:
    #    first_factorisations()
    #else:
    #exit(main())
    N = 1817
    a, b = quadratic_sieve(N)
    print a
    print b

    baseSize = 4
    base = [2, 3, 5, 7]
    dec = np.array([[1, 1, 2, 2], [3, 0, 1, 2], [5, 2, 0, 0], [0, 0, 1, 0], [1, 1, 1, 1]])
    dec_parity = np.array([[1, 1, 0, 0], [1, 0, 1, 0], [1, 0, 0, 0], [0, 0, 1, 0], [1, 1, 1, 1]])
    comb = gauss_jordan(dec_parity, baseSize)
    print comb
    d = np.zeros(baseSize, dtype=np.int)
    for x in comb:
     d[:] += dec[x,:]
    a = 1
    i = 0
    for b in base:
     a *= math.pow(b, d[i]/2)
     i += 1
    g = gcd(a-N, N)
    print g, N/g