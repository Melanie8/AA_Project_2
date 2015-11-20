from sys import stdin, exit
from datetime import datetime
import random
import math
import numpy as np
import copy

print_debug = False
print_time = True
first_factorisations_mode = True
max_prime = 4000
first_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, 2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, 2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741, 2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, 2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, 3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, 3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571, 3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, 3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, 3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, 3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989]
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

def power(x, n):
    n = int(n)
    if n == 0:
        return 1
    if n == 1 or x == 1:
        return x
    y = power(x, n/2)
    if n & 1 == 0:
        return y*y
    else:
        return y*y*x


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


def find_first_primes():
    is_prime = np.ones(max_prime, dtype=bool)
    for i in range(2, max_prime):
        if is_prime[i]:
            first_primes.append(i)
            j = i*i
            while j < max_prime:
                is_prime[j] = False
                j += i


# Compute (n/p)
def legendre(n, p):
    n = n % p
    if n == 0:
        return 0
    ans = 1
    # n =  2^s * m
    s = 0
    m = n
    while m & 1 == 0:
        m /= 2
        s += 1
    if s & 1 == 1 and (p % 8 == 3 or p % 8 == 5):
        ans = -1
    if m == 1:
        return ans
    if m % 4 == 3 and p % 4 == 3:
        ans *= (-1)
    return ans*legendre(p, m)


#For p prime and n with (n/p) = 1, returns r such that r^2 = n mod p
def tonelli_shanks(n, p):
    if p == 2:
        return [1]
    n = n % p
    if p % 4 == 3:
        R = powermod(n, (p+1)/4, p)
        return [R, p-R]
    # p-1 = Q*2^S
    S = 0
    Q = p-1
    while Q & 1 == 0:
        Q /= 2
        S += 1

    z = 2
    while legendre(z, p) != -1:
        z += 1
    c = powermod(z, Q, p)
    R = powermod(n, (Q+1)/2, p)
    t = powermod(n, Q, p)
    M = S
    while True:
        if t == 1:
            return [R, p-R]
        i = 1
        t2i = t*t % p
        while t2i != 1:
            t2i = t2i*t2i % p
            i += 1
        b = powermod(c, int(power(2, M-i-1)), p)
        R = R*b % p
        t = t*b*b % p
        c = b*b % p
        M = i

def square_roots_nul(n, p, q):
    S = tonelli_shanks(n, p)
    R = []
    n = n % q
    for s in S:
        r = s
        while r < q:
            if (r*r)%q == n:
                R.append(r)
            r += p
    return R

def decompose(b, baseSize, base):
    dec = np.zeros(baseSize, dtype=np.long)
    dec_parity = 0
    i = 0
    for f in base:
        if f == -1:
            if b < 0:
                b = -b
                dec[i] = 1
                dec_parity |= 1 << (baseSize-1-i)
        else:
            while b % f == 0:
                dec[i] += 1
                b /= f
            if dec[i] & 1 == 1:
                dec_parity |= 1 << (baseSize-1-i)
        i += 1
    if b > 1:
        return None, None
    else:
        return dec, dec_parity

def gauss_jordan_bitwise(L, nSmooth, baseSize):
    if print_debug:
        print("nSmooth = %d, baseSize = %d" % (nSmooth, baseSize))
    C = [1 << nSmooth-1-i for i in range(nSmooth)]
    indice = [i for i in range(nSmooth)]
    for c in range(baseSize):
        if print_debug:
            print("c = %d" % c)
        # Find pivot
        found = False
        for i in range(c, nSmooth):
            if L[indice[i]] & (1 << (baseSize-c-1)) != 0:
                found = True
                break
        if found:
            temp = indice[i]
            indice[i] = indice[c]
            indice[c] = temp
            if print_debug:
                print("The pivot is in position (%d, %d)" % (i, c))
            # Subtract
            for j in range(c+1, nSmooth):
                if L[indice[j]] & (1 << (baseSize-c-1)) != 0:
                    L[indice[j]] ^= L[indice[c]]
                    C[indice[j]] ^= C[indice[c]]
    return [C[indice[i]] for i in range(nSmooth) if L[indice[i]] == 0]

def gauss_jordan(A, nSmooth, baseSize):
    comb = [[i] for i in range(nSmooth)]
    index = [i for i in range(nSmooth)]
    r = -1
    for j in range(baseSize):
        if print_debug:
            print("A =")
            print A
            print("comb =")
            print comb
        # Find pivot
        for k in range(r+1, nSmooth):
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
                 temp = copy.deepcopy(A[r, :])
                 if print_debug:
                     print temp, A[r, :], A[k, :]
                 A[r, :] = copy.deepcopy(A[k, :])
                 if print_debug:
                     print temp, A[r, :], A[k, :]
                 A[k, :] = copy.deepcopy(temp)
                 if print_debug:
                     print temp, A[r, :], A[k, :]
                 temp = index[k]
                 index[k] = index[r]
                 index[r] = temp
                 if print_debug:
                    print("A =")
                    print A
             for i in range(nSmooth):
                 if i != r:
                     if A[i, j] != 0:
                         for c in comb[index[r]]:
                            comb[index[i]].append(c)
                         for l in range(baseSize):
                             A[i, l] = (A[i, l] - A[r, l]) & 1
                             if print_debug:
                                 print("entry (%d, %d) becomes %d" % (i, l, (A[i, l] - A[r, j]) & 1))
                         if not np.any(A[i, :] % 2):
                             if print_debug:
                                 print "zero line found"
                                 print i
                                 print comb
                                 print comb[index[i]]
                             return comb[index[i]]
    if print_debug:
        print("A =")
        print A
        print("comb =")
        print comb
        print ("index =")
        print index
        print ("returned combination = ")
        print comb[baseSize]
    return comb[baseSize]


def quadratic_sieve(N):
     # Find base
     B = math.ceil(math.pow(math.exp(math.sqrt(math.log(N)*math.log(math.log(N)))), math.sqrt(2)/4))
     if print_debug:
        print("B = %d" % B)
     base = [-1, 2]
     baseSize = 2
     i = 1
     while baseSize < B:
         if print_debug:
            print first_primes[i]
         l = legendre(N, first_primes[i])
         if l == 0:
             return i, N/i
         if l == 1:
             if print_debug:
                print "added"
             base.append(first_primes[i])
             baseSize += 1
         i += 1
     if print_debug:
         print "base = ",
         print base

     # Find >baseSize smooth numbers in the interval  [sqrt(n)-M, sqrt(n)+M] with M = B^3
     M = int(power(B, 3))
     X = max(1, int(math.ceil(math.sqrt(N)-M)))


     #method = 1
     #method = 2
     method = 2 # do both

     if method != 1:
         T = np.zeros(2*M, dtype=np.long)
         smooth = []
         nSmooth = 0
         log_base = [math.ceil(math.log(base[i])) for i in range(1, baseSize)]
         for p in range(1, baseSize):
             q = base[p]
             while q <= 100000: # Now also do the powers of p
                 if q == base[p]:
                     R = tonelli_shanks(N, q)
                 else:
                     R = square_roots_nul(N, base[p], q) # If q = p^j with j >= 2, our technique for finding the square roots is a bit silly
                     if len(R) == 0: # If there are not square roots mod q = p^j then there will be no square roots mod p^(j+1)
                         break
                 print("q = %d" % q)
                 if print_debug:
                     print("q = %d, R = " % q),
                     print R
                 m = X % q
                 for r in R:
                     if r-m >= 0:
                         i = r-m
                     else:
                         i = r-m+q
                     i = int(i)
                     while i < 2*M:
                         T[i] += log_base[p-1]
                         i += q
                 q *= base[p]
         #target = math.log(N)/2+math.log(M)
         #thresh = 3/4*target
         E = [math.exp(i) for i in range(1, int(math.floor(math.log(M))) + 2)] # Compute e^1, e^2, e^3, ...
         target2 = [math.log(N)/2+i+1 for i in range(1, int(math.floor(math.log(M))) + 2)] # Compute the different targets (will depend on i)
         thresh2 = [3*i/4 for i in target2]
         current = 0
         # Guys >= sqrt(n)
         for j in range(0,M):
             if E[current+1] <= j:
                 current = current+1
             i = M+j
             if T[i] >= thresh2[current]:
                 smooth.append(X+i)
                 nSmooth += 1
         current = 0
         # Guys < sqrt(n)
         for j in range(1,M):
             if E[current+1] <= j:
                 current = current+1
             i = M-j
             if T[i] >= thresh2[current]:
                 smooth.append(X+i)
                 nSmooth += 1

     if method == 3:
         T2 = T

     if method != 2:
         T = np.zeros(2*M, dtype=np.long)
         #smooth = []
         #nSmooth = 0
         for i in range(2*M):
             T[i] = abs((X+i)*(X+i)-N)
             #if T[i] == 1:
             #    smooth.append(X+i)
             #    nSmooth += 1
             #    if print_debug:
             #        print("!!!! T2[%d] = %d, target = %d" % (i, T2[i], target))
             if T[i] == 0:
                 return X+i, X+i
         if print_debug:
             print "T = ",
             print T
         for p in range(1, baseSize):
             R = tonelli_shanks(N, base[p])
             print("p = %d" % base[p])
             if print_debug:
                 print("p = %d, R = " % base[p]),
                 print R
             m = X % base[p]
             for r in R:
                 if r-m >= 0:
                     i = r-m
                 else:
                     i = r-m+base[p]
                 i = int(i)
                 while i < 2*M:
                     if T[i] % base[p] != 0:
                         print("WTF!!! p = %d, r = %d, i = %d, T[i] = %d" % (base[p], r, i, T[i]))
                     T[i] /= base[p]
                     while T[i] % base[p] == 0:
                         T[i] /= base[p]
                     #if T[i] == 1:
                     #    smooth.append(X+i)
                     #    nSmooth += 1
                     #    if print_debug and method == 3:
                     #       print("!!!! T2[%d] = %d, target = %d" % (i, T2[i], target))
                     i += base[p]
     
     '''
     if method == 3:
         falsesmooth = 0
         smoothmissed = 0
         realsmooth = 0
         for i in range(2*M):
             if T[i] == 1 and T2[i] >= 75:
                 realsmooth += 1
             elif T[i] == 1:
                 smoothmissed += 1
             elif T2[i] >= 75:
                 falsesmooth += 1          
         print("REAL SMOOTH : %d", realsmooth)
         print("FALSE SMOOTH : %d", falsesmooth)
         print("MISSED SMOOTH : %d", smoothmissed)
     '''
         
     print("SMOOTHS FOUND BEFORE VERIF : %d, BASE SIZE : %d" % (nSmooth, baseSize))

     # Decompose them in the base
     dec = []
     dec_parity = []
     real_smooth = []
     for i in range(nSmooth):
         if print_debug:
            print("smooth[i] = %d" % smooth[i])
         a, b = decompose(smooth[i]*smooth[i]-N, baseSize, base)
         if a is None:
             nSmooth -= 1
             if print_debug:
                 print "not a smooth number"
         else:
             if print_debug:
                 print a
                 print "{0:b}".format(b)
             dec.append(a)
             dec_parity.append(b)
             real_smooth.append(smooth[i])
     print("SMOOTHS FOUND AFTER VERIF : %d" % nSmooth)
     
     if nSmooth <= baseSize:
         print("Not enough smooth numbers :(")

     # Find even combination
     combs = gauss_jordan_bitwise(dec_parity, nSmooth, baseSize)
     if print_debug:
         print "combs = ",
         print combs

     # Find non trivial factors
     for comb in combs:
         a = 1
         b = 1
         tot = np.zeros(baseSize, dtype=np.int)
         if print_debug:
             print "comb = ",
             print "{0:b}".format(comb)
         x = nSmooth-1
         while comb != 0:
             if comb & 1 == 1:
                 a = (a*real_smooth[x])%N
                 for i in range(baseSize):
                     if dec[x][i] != 0:
                        tot[i] += dec[x][i]
             comb = comb >> 1
             x -= 1
         for i in range(baseSize):
             b = (b * powermod(base[i], tot[i]/2, N)) % N
         if print_debug:
             print a, b
         g = gcd(a+b, N)
         if g != 1 and g != N:
             return g, N/g
     return -1, -1


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
    #find_first_primes()
    time0 = datetime.now()
    N = 138762438070143122474623
    #138762438070143122474623 #130750493899956009090497 #4153755852461398439 #32193886486049401L #253442228047 #152061056489 #152020131149 #7373 #5535333969127 #3691153417681 #10379
    a, b = quadratic_sieve(N)
    print a
    print b
    elapsed_time_function("main", time0)
    '''a, b = decompose(84, 4,[2, 3, 5, 7])
    print a
    print b
    print "{0:b}".format(b)'''


