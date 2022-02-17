from math import sqrt, isqrt, prod, log, exp, floor, gcd
from utils import legendre_symbol, valuation
from prime_sieves import sieve_of_Eratosthenes
from random import randint
from Crypto.Util.number import *
import numpy as np
from tqdm import tqdm

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
        a = round((x + y) / z)

        Ak = (e2 + x*f2) % n
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
    Xis = [0]*(m+2)
    x_set = set()
    matrix = np.zeros((m+2, m+1), dtype=np.int8)
    count = 0
    with tqdm(total=m+2) as pbar:
        for i, j in enumerate(quadratic_convergents_numerator(n)):
            if count == m+2:
                break
            t = (j*j) % n
            #  P, Q = j
            #  t = P**2 - n*Q**2
            if is_prod_of_basis(t, k):
                #  print("True")
                if j in x_set:
                    continue
                Xis[count] = j
                x_set.add(j)
                for kk, l in enumerate(basis):
                    matrix[count, kk+1] = valuation(t, l) % 2
                count+=1
                pbar.update(1)
            elif is_prod_of_basis(-t % n, k):
                #  print("True")
                if j in x_set:
                    continue
                Xis[count] = j
                x_set.add(j)
                t = -t % n
                matrix[count, 0] = 1
                for kk, l in enumerate(basis):
                    matrix[count, kk+1] = valuation(t, l) % 2
                count+=1
                pbar.update(1)
    B = np.zeros(m+1)
    return matrix.transpose()
    res = np.linalg.lstsq(matrix.transpose(), B, rcond=None)[0] % 2
    X = 0
    for i, j in zip(Xis, res):
        X += i*j
    print(X)
    

def test(nb):
    n = getPrime(nb)*getPrime(nb)
    factor_cfrac(n)

test(64)
