from math import sqrt, isqrt, prod, log, exp, floor, gcd
from utils import legendre_symbol, valuation, is_square
from prime_sieves import sieve_of_Eratosthenes
import numpy as np
from tqdm import tqdm
from gauss_elim import find_perfect_square, try_combination
from typing import Generator, Optional

def generate_factor_basis(n:int) -> list[int]:
    """
    Generates a list of primes to be used as basis in factoring using CFRAC
    The smoothness bound B is computed using an optimised value.

    :param n: The int being factored
    :return : List of prime forming the base
    """
    temp = sqrt(log(n) * (log(log(n)))) / 2
    B = floor(exp(temp))
    primes = sieve_of_Eratosthenes(B)
    return [2] + [p for p in primes[1:] if legendre_symbol(n, p)==1]


def is_prod_of_basis(n:int, k:int) -> bool:
    """
    Checks if a given number is smooth over a factor base

    :param n: The number to be checked
    :param k: Product of primes in the base
    :return : Boolean
    """
    g = gcd(n, k)
    while g>1:
        while n%g == 0:
            n//=g
        if n == 1:
            return True
        g = gcd(n, k)
    return False


def quadratic_convergents_numerator(n:int) -> Generator[int, int, None]:
    """
    Computes Continued fraction convergents of sqrt(n)

    reference: https://trizenx.blogspot.com/2018/10/continued-fraction-factorization-method.html

    :param n: target of square root continued fraction approximations
    :return: A generator giving only the numerator of the convergents, uncomment code for denominator
    """
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

def generate_relations(n:int, primes:list[int], mul:int) -> tuple[list[int], np.array]:
    """
    Finds all the relations of x^2 = t mod n, where t is smooth over given basis using CFRAC
    This step can be slow for larger numbers so tqdm is included as a progress bar.

    :param n: The number being factored
    :param primes: The factor basis, a list of primes
    :param mul: Optional multiplier used for convergents incase using just n fails
    :return : List of all x, matrix containing exponents vectors of t over the basis
    """
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

def factor_cfrac(n:int) -> Optional[list[tuple[int, int]]]:
    """
    Returns list of prime factors of the given semiprime using the CFRAC method
    
    :param n: semiprime being factored
    :return: list of factors
    """
    print("[+] Given n:",n)

    if is_square(n):
        print("[+] n is a perfect square")
        facs = [(isqrt(n), 2)]
        print("\n[+] Factors found:", facs)
        return facs

    print("[+] Generating factor basis")
    primes = generate_factor_basis(n)
    for mul in [1, 2, 3, 5, 7, 11, 13]:
        print("\n[+] Generating realtions with multiplier", mul)
        Xis, matrix = generate_relations(n, primes, mul)
        matrix_backup = matrix.tolist()
        for combination in find_perfect_square(matrix):
            facs = try_combination(Xis, primes, matrix_backup, n, combination)
            if facs:
                print("\n[+] Factors found:", facs)
                return facs
        else:
            print("[+] No factors found")
    return None
