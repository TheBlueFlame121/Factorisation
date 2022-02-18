import numpy as np
from numba import njit
from math import gcd
from typing import Generator, Optional

"""
Helper file containing all the required functions for Gaussian Elimination step in
Continued Fraction method, Quadratic Sieve and Number Field Sieve.

After all relations of the form x^2 = t mod N, where t is smooth have been found, 
find_perfect_square will give combinations to form X^2 = Y^2 mod N and
try_combination will try test that combination
"""


@njit
def fast_gauss(matrix:np.array) -> tuple[list[bool], np.array]:
    """
    Applies Fast Gaussian elimination to a binary (GF(2)) matrix.
    Uses numba for speedup and hence is implemented in a more basic manner

    reference: https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf
                
    :param matrix: A numpy matrix containing the exponents vectors of all the the
    :return : A list of marked rows and a gaussian eliminated matrix
    """
    matrix%=2
    m, n = matrix.shape
    mark = [False]*m

    for j, column in enumerate(matrix.T):
        i=0
        while i < n  and matrix[i, j] != 1:
           i+=1
        if i==n:
            continue
        pivot = i
        mark[pivot] = True

        for k, col in enumerate(matrix.T):
            if k == j:
                continue
            if col[pivot] == 1:
                matrix[:, k] += matrix[:, j]
                matrix[:, k] %= 2

    return mark, matrix

def find_perfect_square(matrix:np.array) -> Generator[list[int], np.array, None]:
    """
    Uses fast gauss to find and return valid combinations of x^2 = t mod N such that
    their product is of the form X^2 = Y^2 mod N.

    reference: https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf

    :param matrix: List of exponent vectors of all t
    :return : Valid combination as a list of row indices
    """
    print("[+] Running gauss elimination")
    mark, m = fast_gauss(matrix)
    m = m.tolist()

    independent = [(r.index(1), i,r) for i, r in enumerate(m) if mark[i]]
    dependent = [(i,r) for i, r in enumerate(m) if not mark[i]]
    
    print("[+] Trying combinations")
    for index, row in dependent:
        res = [index]

        check = 0
        if all(i==0 for i in row):
            continue

        for index, element in enumerate(row):
            if element == 1:
                for one, i, r in independent:
                    if one == index:
                        res.append(i)
                        break
                else:
                    check = 1
                    break
        if check == 1:
            continue
        yield res            

def try_combination(Xis:list[int], primes:list[int], matrix:list[list[int]], n:int, comb:list[int]) -> Optional[list[int]]:
    """
    Given a valid combination of rows, computes X, Y and tries to factor N

    :param Xis: list of all x from relations
    :param primes: The factor basis
    :praram matrix: All the expoenent vectors of t
    :param n: The number to be factored
    :param combination: A valid combination of rows such that the product of ts is a square
    :return : Factors of n if found, else None
    """
    x = 1
    m = len(primes)
    primes = [-1] + primes
    
    exponents = np.zeros(m+1, dtype=np.int32)

    for i in comb:
        exponents += np.array(matrix[i])
        x *= Xis[i]
        x %= n

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
