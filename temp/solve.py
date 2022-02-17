from math import sqrt, isqrt, prod, log, exp, floor, gcd
from utils import legendre_symbol, valuation
from prime_sieves import sieve_of_Eratosthenes
from random import randint
from Crypto.Util.number import *
import numpy as np
from tqdm import tqdm
from gauss_elim import find_perfect_square

def generate_factor_basis(n):
    temp = sqrt(log(n) * (log(log(n)))) / 2
    B = floor(exp(temp))
    primes = sieve_of_Eratosthenes(B)
    return [2] + [p for p in primes[1:] if legendre_symbol(n, p)==1]


def is_prod_of_basis(n, k):
    g = gcd(n, k)
    while g>1:
        while n%g == 0:
            n//=g
        if n == 1:
            return True
        g = gcd(n, k)
    return False


def quadratic_convergents_numerator(n):
    x = y =  isqrt(n)
    z = 1
    a = 2*x

    e1, e2 = 1, 0
    f1, f2 = 0, 1

    while True:
        y = a*z - y
        z = (n - y*y) // z
        a = round((x + y) / z)

        Ak = (e2 + x*f2) % n
        #  Bk = f2

        yield Ak#, Bk

        f1, f2 = f2 % n, (a*f2 + f1) % n
        e1, e2 = e2 % n, (a*e2 + e1) % n

        if z == 1:
            break

def generate_base(n, primes, mul):
    m = len(primes)
    k = prod(primes)

    Xis = [0]*(m+2)
    x_set = set()
    matrix = np.zeros((m+2, m+1), dtype=np.int32)
    count = 0
    pbar = tqdm(total=m+2, unit=" reln", desc="Progress ")
    for i, j in enumerate(quadratic_convergents_numerator(mul*n)):
        if count == m+2:
            break
        t = (j*j) % n
        if is_prod_of_basis(t, k):
            if j in x_set:
                continue
            Xis[count] = j
            x_set.add(j)
            for kk, l in enumerate(primes):
                matrix[count, kk+1] = valuation(t, l)
            count+=1
            pbar.update(1)
        elif is_prod_of_basis(-t % n, k):
            if j in x_set:
                continue
            Xis[count] = j
            x_set.add(j)
            t = -t % n
            matrix[count, 0] = 1
            for kk, l in enumerate(primes):
                matrix[count, kk+1] = valuation(t, l)
            count+=1
            pbar.update(1)
    return Xis, matrix

def try_combination(Xis, primes, matrix, n, comb):
    x = 1
    for i in comb:
        x *= Xis[i]
        x %= n

    m = len(primes)
    primes = [-1] + primes
    
    exponents = np.zeros(m+1, dtype=np.int32)

    for i in comb:
        exponents += np.array(matrix[i])

    exponents//=2
    
    y = 1
    for i, j in zip(primes, exponents):
        y *= pow(i, int(j), n)
        y %= n

    assert pow(x, 2, n) == pow(y, 2, n), print(combination)
    if x == y or x == -y%n:
        return None

    factor = gcd(x-y, n)
    if factor == 1:
        return None
    return [(factor, 1), (n//factor, 1)]

def factor_cfrac(n):
    print("[+] Generating factor basis")
    primes = generate_factor_basis(n)
    print("\n[+] Generating Base matrix")
    Xis, matrix = generate_base(n, primes)
    count = 0
    for combination in find_perfect_square(matrix):
        print(count)
        count+=1
        facs = try_combination(Xis, primes, matrix, n, combination)
        if facs:
            return facs
    
def test(nb):
    n = getPrime(nb)*getPrime(nb)
    print(factor_cfrac(n))

#  test(64)
n = getPrime(64)*getPrime(64)
# Remove factors of 2 and check if perfect square as well
print(n)
primes = generate_factor_basis(n)
#  print(primes)
for mul in [1, 3, 5, 7, 11, 13]:
    Xis, matrix = generate_base(n, primes, mul)
    #  for i, j in zip(Xis, matrix):
        #  print(i, j)
    matrix_backup = matrix.tolist()
    for combination in find_perfect_square(matrix):
        facs = try_combination(Xis, primes, matrix_backup, n, combination)
        if facs:
            print(facs)
            exit()
