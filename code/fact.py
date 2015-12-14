from sys import stdin, exit
from datetime import datetime
import random
import math
import numpy as np
import copy

first_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, 2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, 2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741, 2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, 2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, 3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, 3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571, 3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, 3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, 3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, 3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989]

print_debug = False
print_time = False

# Print the elapsed time between time0 and now
def elapsed_time_function(function_name, time0):
    time1 = datetime.now()
    elapsed_time = time1 - time0
    s = 'Elapsed time %s: %d seconds and %d milliseconds' % (function_name, elapsed_time.seconds, elapsed_time.microseconds/1000)
    if print_time:
        print(s)
    return [elapsed_time.seconds, elapsed_time.microseconds/1000]


# Compute the GCD of a and b
def gcd(a, b):
    while b != 0:
        tmp = b
        b = a % b
        a = tmp
    return a


# Compute x^n % m efficiently
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


# Compute x^n efficiently
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


# Miller Rabin test : is N a prime number?
def miller_rabin(N):
    if print_debug:
        print("Miller Rabin for N = %d" % N)
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
    for i in range(20): # 20 tests : Probability of false prime = 1/4^20
        prob_prime = False
        a = random.randrange(1, N)
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

######################################
#           Pollard Rho              #
######################################
# Pollard Rho algorithm to factorize N
def pollard_rho(N, param):
    i=1
    if print_debug:
        print("Pollard rho for N = %d" % N)
    x = 2
    saved_x = 2
    y = 2
    saved_y = 2
    r = 1
    d = 1
    e = 0
    while d == 1:
        i+=1
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
            i+=1
            x = (x*x + param) % N
            temp = (y*y + param) % N
            y = (temp*temp + param) % N
            if x > y:
                d = gcd(x-y, N)
            else:
                d = gcd(y-x, N)
        if d == N:
            return pollard_rho(N, param + 1)
    return d, N/d,i


# Brent version of Pollard Rho algorithm to factorize N
def brent(N, c):
    if print_debug:
        print("Brent for N = %d" % N)
    if N % 2 == 0:
        return 2, N/2
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
            y = (y*y+c) % N
        while True:
            ys = y
            for i in range(min(m, r-k)):
                q = q*abs(x-y) % N
                y = (y*y+c) % N
            g = gcd(q, N)
            k += m
            if k >= r or g > 1:
                break
        r *= 2
        if g > 1:
            break
    
    if g == N:
        while True:
            g = gcd(abs(x-ys),N)
            ys = (ys*ys+c) % N
            if g > 1:
                break
        if g == N:
            return brent(N, c+1)
    return g, N/g


######################################
#               Fermat               #
######################################
# Fermat 1 to factorize N
def fermat1(N):
    if print_debug:
        print("Fermat 1 for N = %d" % N)
    for a in range(int(math.ceil(math.sqrt(N))), (N+9)/6):
        b = math.sqrt(a*a-N)
        if math.floor(b) == b:
            b = int(b)
            return a-b, N/(a-b)


# Fermat 2 to factorize N
def fermat2(N):
    if print_debug:
        print("Fermat 2 for N = %d" % N)
    a = int(math.ceil(math.sqrt(N)))
    b = a*a-N
    sqrtb = math.sqrt(b)
    while math.ceil(sqrtb) != sqrtb:
        b += 2*a+1
        a += 1
        sqrtb = math.sqrt(b)
    return a-int(sqrtb), a+int(sqrtb)


######################################
#          Quadratic Sieve           #
######################################
# Compute the Legendre symbol (n/p)
def legendre(n, p):
    if print_debug:
        print("Compute (%d/%d)" % (n, p))
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


# For p prime and n with (n/p) = 1, returns all r such that r^2 = n mod p
def tonelli_shanks(n, p):
    if print_debug:
        print("Tonelli Shanks for finding square roots of %d mod %d" % (n, p))
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


# For q = p^k with k >= 2 and n with (n/p) = 1, returns all r such that r^2 = n mod q
def square_roots_bad(n, p, q):
    if print_debug:
        print("Finding square roots of %d mod %d" % (n, q))
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


# Decompose number b with the primes of the base, if this is possible
def decompose(b, baseSize, base):
    if print_debug:
        print("Decompose %d" % b)
    if b == 0:
        return -1, -1
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
    # b cannot be factorized in the base
    if b > 1:
        return None, None
    else:
        return dec, dec_parity


# Gauss Jordan algorithm
def gauss_jordan_bitwise(L, nSmooth, baseSize):
    if print_debug:
        print "Gauss-Jordan"
    C = [1 << nSmooth-1-i for i in range(nSmooth)]
    indice = [i for i in range(nSmooth)]
    for c in range(baseSize):
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
            # Subtract
            for j in range(c+1, nSmooth):
                if L[indice[j]] & (1 << (baseSize-c-1)) != 0:
                    L[indice[j]] ^= L[indice[c]]
                    C[indice[j]] ^= C[indice[c]]
    return [C[indice[i]] for i in range(nSmooth) if L[indice[i]] == 0]


# Quadratic Sieve algorithm to factorize N
def quadratic_sieve(N):
     if print_debug:
        print("Quadratic Sieve to factorize %d" % N)
     # Find base
     B = math.ceil(math.pow(math.exp(math.sqrt(math.log(N)*math.log(math.log(N)))), math.sqrt(2)/4))
     base = [-1, 2]
     baseSize = 2
     i = 1
     while baseSize < B:
         l = legendre(N, first_primes[i])
         if l == 0:
             return i, N/i
         if l == 1:
             base.append(first_primes[i])
             baseSize += 1
         i += 1

     M = int(power(B*1, 3))/3
     # Find >baseSize smooth numbers in the interval  [sqrt(n)-M, sqrt(n)+M]
     X = max(1, int(math.ceil(math.sqrt(N)-M)))

     T = np.zeros(2*M, dtype=np.int)
     smooth = []
     nSmooth = 0
     for p in range(1, baseSize):
         log = int(math.ceil(math.log(base[p])))
         q = base[p]
         while q <= 100000: # Now also do the powers of p
             if q == base[p]:
                 R = tonelli_shanks(N, q)
             else:
                 R = square_roots_bad(N, base[p], q) # If q = p^j with j >= 2, our technique for finding the square roots is a bit silly
                 if len(R) == 0: # If there are not square roots mod q = p^j then there will be no square roots mod p^(j+1)
                     break
             if print_debug:
                 print("q = %d" % q)
             m = X % q
             for r in R:
                 if r-m >= 0:
                     i = r-m
                 else:
                     i = r-m+q
                 i = int(i)
                 while i < 2*M:
                     T[i] += log
                     i += q
             q *= base[p]
     
     E = [math.exp(i) for i in range(1, int(math.floor(math.log(M))) + 5)] # Compute e^1, e^2, e^3, ...
     target2 = [math.log(N)/2+i+1 for i in range(1, int(math.floor(math.log(M))) + 5)] # Compute the different targets (will depend on i)
     thresh2 = [int(3*i/4) for i in target2]
     current = 0
     # Guys >= sqrt(n)
     j = 0
     while j < M:
         if E[current+1] <= j:
             current = current+1
         i = M+j
         if T[i] >= thresh2[current]:
             smooth.append(X+i)
             nSmooth += 1
         j += 1
     current = 0
     # Guys < sqrt(n)
     j = 1
     while j < M:
         if E[current+1] <= j:
             current = current+1
         i = M-j
         if T[i] >= thresh2[current]:
             smooth.append(X+i)
             nSmooth += 1
         j += 1
     if print_debug:
         print("Potential smooth numbers found : %d (Base size : %d)" % (nSmooth, baseSize))
     
     # Decompose in the base
     dec = []
     dec_parity = []
     real_smooth = []
     for i in range(nSmooth):
         a, b = decompose(smooth[i]*smooth[i]-N, baseSize, base)
         if a is None:
             nSmooth -= 1
         elif b == -1:
             return smooth[i], smooth[i] # Particular case of a square
         else:
             dec.append(a)
             dec_parity.append(b)
             real_smooth.append(smooth[i])
     if print_debug:
         print("Real smooth numbers found : %d" % nSmooth)
     
     if print_debug:
         if nSmooth <= baseSize+1:
             print("Not enough smooth numbers :(")

     # Find even combinations
     combs = gauss_jordan_bitwise(dec_parity, nSmooth, baseSize)
     if print_debug:
         print "combs = ",
         print combs

     # Find non trivial factors
     for comb in combs:
         a = 1
         b = 1
         tot = np.zeros(baseSize, dtype=np.int)
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
         g = gcd(a+b, N)
         if g != 1 and g != N:
             return g, N/g
     return 1, N


def xGCD(a, b):
    if b == 0:
       x = 1
       y = 0
       return x

    x = xGCD(b, a % b)
    x = y1
    y = x1 - (a / b) * y1
    return x


# Factorize N
def factorize(N):
    if print_debug:
        print("Factorize %d" % N)
    factors = []
    tofactor = []
    while N & 1 == 0:
        N /= 2
        factors.append(2)
    if miller_rabin(N):
        factors.append(N)
    else:
        tofactor.append(N)
        while(tofactor):
            newtofactor = []
            for f in tofactor:
                p0, p1 = quadratic_sieve(f)
                if p0 == 1 or p1 == 1:
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
    time0 = datetime.now()
    userinput = int(stdin.readline())
    while userinput!=0:
         if print_debug:
            print('\nInput = %s\n' % userinput)
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


# Factorize all the small numbers (not used anymore)
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
    '''time0 = datetime.now()
    N = 152020131149
    print N, math.pow(N, 0.25)
    a,b,c=pollard_rho(N,1)
    print a, math.sqrt(a)
    print b
    print c
    d,e= elapsed_time_function("main", time0)
    print d
    print e
    print c/(d+0.001*e)'''

    exit(main())
    
    # Some statistics about quadratic sieve:
    #23265811160790642537889302884495500847 --> 2040 sec (B habituel -> 1015, M = B^3/3 -> ?, SMOOTHS : 300948 -> 4343 / 1015)
    #4458135213943293827010603080147 --> 228 sec (B habituel -> 459, M = B^3/3 -> 32234193, SMOOTHS : 48756 -> 1067 / 459)
    #83881350134177350355462809427 --> 117 sec (B habituel -> 370, M = B^3/3 -> 16884333, SMOOTHS : 28035 -> 558 / 370)
    #46047434098244260962615985003 --> 162 sec (B habituel -> 358, M = B^3/3 -> 15294237, SMOOTHS : 68080 -> 1251 / 358)
    #460474340983462127457456719 --> 47 sec (B habituel -> 277, M = B^3/3 -> 7084644, SMOOTHS : 14126 -> 363 / 277)
    #138762438070143122474623 --> 13 sec (B habituel -> 172, M = B^3/3 -> 1696149, SMOOTHS : 8912 -> 211 / 172)
    #130750493899956009090497 --> 15 sec (B habituel -> 172, M = B^3/3 -> 1696149, SMOOTHS : 7634 -> 243 / 172)
    #4153755852461398439 --> 2 sec (B habituel -> 89, M = B^3/3 -> 234989, SMOOTHS : 6172 -> 376 / 89)
    #32193886486049401 --> 1 sec (B habituel -> 64, M = B^3/3 -> 87381, SMOOTHS : 18599 -> ? / 64)
    #5535333969127 --> 0.146 sec (B habituel -> 34, M = B^3/3 -> 13101, SMOOTHS : 730 -> 44 / 34)
    #3691153417681 --> 0.317 sec (B habituel -> 33, M = B^3/3 -> 11979, SMOOTHS : 2534 -> 221 / 33)
    #253442228047 --> 0.089 sec (B habituel -> 27, M = B^3/3 -> 6561, SMOOTHS : 598 -> 57 / 27)
    #152061056489 --> 0.108 sec (B habituel -> 26, M = B^3/3 -> 5858, SMOOTHS : 777 -> 81 / 26)
    #152020131149 --> 0.065 sec (B habituel -> 26, M = B^3/3 -> 5858, SMOOTHS : 357 -> 31 / 26)
    #10379 --> 0.009 sec (B habituel -> 5, M = B^3/3 -> 41, SMOOTHS : 11 -> 4 / 5)
