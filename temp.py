from math import sqrt, isqrt, prod, log, exp, floor, gcd
from utils import legendre_symbol, valuation
from prime_sieves import sieve_of_Eratosthenes
from random import randint
from Crypto.Util.number import *

def generate_factor_basis(n):
    temp = sqrt(log(n) * (log(log(n)))) / 2
    B = floor(exp(temp))
    primes = sieve_of_Eratosthenes(B)
    return [p for p in primes if legendre_symbol(n, p)==1]


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
        a = (x + y) // z

        Ak = (e2 + x*f2)
        #  Bk = f2

        yield Ak#, Bk

        f1, f2 = f2 % n, (a*f2 + f1) % n
        e1, e2 = e2 % n, (a*e2 + e1) % n

        if z == 1:
            break

def factor_cfrac(n):
    basis = generate_factor_basis(n)
    print("N modulo 4 = ", n%4)
    m = len(basis)
    print(m)
    k = prod(basis)
    Xis = []
    matrix = []
    for i, j in enumerate(quadratic_convergents_numerator(n)):
        if i == 10000:
            break
        t = (j*j) % n
        #  P, Q = j
        #  t = P**2 - n*Q**2
        if is_prod_of_basis(t, k):
            print("True")
            Xis.append(j)
            powers = [0]*(m+1)
            for k, l in enumerate(basis):
                powers[k+1] = valuation(t, l) % 2
            matrix.append(powers)
        elif is_prod_of_basis(-t % n, k):
            print("True")
            Xis.append(j)
            t = -t % n
            powers = [1] + [0]*m
            for k, l in enumerate(basis):
                powers[k+1] = valuation(t, l) % 2
            matrix.append(powers)
    print(len(Xis))
