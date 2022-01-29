from utils import valuation
from math import isqrt

def sieve_of_Eratosthenes(n:int) -> list[int]:
    nums = [True for i in range(n+1)]
    primes = []
    p = 2
    while (p**2 <= n):
        if (nums[p] == True):
            for i in range(p**2, n+1, p):
                nums[i] = False
        p += 1
    for p in range(2, n+1):
        if nums[p]:
            primes.append(p)
    return primes

def naive_Sundaram(n: int) -> list[int]:
    primes = [2]
    k = (n - 2) // 2
    integers_list = [True] * (k + 1)
    for i in range(1, k + 1):
        j = i
        while i + j + 2*i*j <= k:
            integers_list[i + j + 2 * i * j] = False
            j += 1
    for i in range(1, k + 1):
        if integers_list[i]:
            primes.append(2*i + 1)
    return primes

def sieve_of_Sundaram(n: int) -> list[int]:
    if n < 3:
        if n < 2: return []
        else: return [2]    
    k = (n - 3) // 2 + 1
    integers_list = [True] * k
    primes = [2]
    for i in range(0, (isqrt(n) - 3) // 2 + 1):
            p = 2 * i + 3
            s = (p * p - 3) // 2 
            for j in range(s, k, p):
                integers_list[j] = False
    for i in range(0, k):
        if integers_list[i]:
            primes.append(2*i + 3)
    return primes

def factor_erat(n:int) -> list[tuple]:
    primes = sieve_of_Eratosthenes(n)
    fac = []
    for p in primes:
        if n%p == 0:
            fac.append((p, valuation(n, p)))
    return fac

def factor_sund(n:int) -> list[tuple]:
    primes = sieve_of_Sundaram(n)
    fac = []
    for p in primes:
        if n%p == 0:
            fac.append((p, valuation(n, p)))
    return fac
