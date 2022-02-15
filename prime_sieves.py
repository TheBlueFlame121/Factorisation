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

# Utility functions for the bitsieve
def getmask(n, i):
    x, y = divmod(n, i)
    return x, 1<<(i-y-1)

def int_to_primes(n, i, b):
    res = []
    offset = i*b
    for i, j in enumerate(bin(n)[2:].zfill(b)):
        if j=='1':
            res.append(offset+i)
    return res

# Eratosthenes seive implemented using bits instead of booleans
# Slightly slower but much more efficient memory wise
# To account for slower speed I tried ctypes int32 but that was slower still
def bitsieve(n:int, b:int=30) -> list[int]:
    x, y = divmod(n, b)
    primes = [2**b - 1]*(x+1)
    primes[-1] ^= (1<<(b-y-1)) - 1
    primes[0] ^= ((1<<(b-2)) | (1<<(b-1)))
    for p in range(2, isqrt(n)+1):
        x, y = getmask(p, b)
        if (primes[x] & y):
            for i in range(p**2, n+1, p):
                x, y = getmask(i, b)
                primes[x] = primes[x] & ~y
    res = []
    for i, j in enumerate(primes):
        res += int_to_primes(j, i, b)
    return res

# Dumbed down version of the Sundaram Sieve
# https://en.wikipedia.org/wiki/Sieve_of_Sundaram
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

# Actual Sundaram Sieve
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

# Algorithm from Wikipedia
# https://en.wikipedia.org/wiki/Sieve_of_Atkin
def sieve_of_Atkin(limit: int) -> list[int]:
    # Integers mod 60 with wheel 2/3/5
    # It's faster to just store them directly
    s = [1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 54, 59]

    primes = [False]*limit
    for w in range(0, limit//60):
        for x in s:
            primes[60*w + x] = False

    x = 1
    while x*x < limit:
        y = 1
        while y*y <= limit:
            n = 4*x*x + y*y
            if n<limit and n%60 in [1, 13, 17, 29, 37, 41, 49, 53]:
                primes[n]= primes[n]^True
            y+=2
        x+=1

    x = 1
    while x*x <= limit:
        y = 2
        while y*y < limit:
            n = 3*x*x + y*y
            if n<limit and n%60 in [7, 19, 31, 43]:
                primes[n] = primes[n]^True
            y+=2
        x+=2

    x = 1
    while x*x <= limit:
        y = x-1
        while y > 0:
            n = 3*x*x - y*y
            if n>0 and n<limit and n%60 in [11, 23, 47, 59]:
                primes[n] = primes[n]^True
            y-=2
        x+=1

    # This is what's mentioned of wikipedia but it lets some composites through
    #  for n in [60*w + x for x in s for w in range(limit)]:
    #      if n<7:
    #          continue
    #      if n*n > limit:
    #          break
    #      if primes[n]:
    #          for c in [n*n*(60*w + x) for x in s for w in range(limit)]:
    #              if c>limit:
    #                  break
    #              primes[c]=False

    # Alternative implementation but less effective
    r = 5
    while r * r < limit:
        if primes[r]:
            for i in range(r * r, limit, r * r):
                primes[i] = False
        r += 1

    res = [2, 3, 5]
    n = 7
    while n < limit:
        if primes[n]:
            res.append(n)
        n+=1

    return res

def factor_sieve(n:int, sieve) -> list[tuple[int, int]]:
    primes = sieve(n)
    fac = []
    for p in primes:
        if n%p == 0:
            fac.append((p, valuation(n, p)))
    return fac
