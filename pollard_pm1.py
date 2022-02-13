from math import gcd, floor, log
from prime_sieves import sieve_of_Eratosthenes
from sympy import prevprime
from utils import valuation
from random import randint

# Followed the Algorithm from https://en.wikipedia.org/wiki/Pollard%27s_p_%E2%88%92_1_algorithm
# Assumes input is semiprime
# It's a probabilistic algorithm so it doesn't always succeed
# Commented out lines are the naive approaches instead of the optimised ones
def factor(n:int, B1:int, B2:int) -> list[tuple]:
    primes = sieve_of_Eratosthenes(B2)
    ind = primes.index(prevprime(B1)) + 1
    small_primes, large_primes = primes[:ind], primes[ind:]

    #  M=1
    #  for q in small_primes:
    #      M*=pow(q, floor(log(B1, q)))

    for _ in range(4):
        a = randint(2, n)

        if gcd(a, n) != 1:
            g = gcd(a, n)
            v = valuation(n, g)
            return [(g, v), (n//(g**v), 1)]      
        
        #  H = pow(a, M, n)
        H = a
        for q in small_primes:
            for i in range(floor(log(B1, q))):
                H = pow(H, q, n)
        
        g = gcd(H-1, n)

        if g in [1, n]:
            #  Q = 1
            #  for q in large_primes:
            #      Q*= pow(H, q, n) - 1
            #  g = gcd(Q, n)
            
            diff = [j-i for i, j in zip(large_primes, large_primes[1:])]
            table = dict()
            temp = pow(H, 2, n)
            for i in range(2, max(diff)+2, 2):
                table[i] = temp
                temp *= pow(H, 2, n)
                temp %= n

            temp2 = pow(H, large_primes[0], n)
            Q = temp2 - 1
            for q1, q2 in zip(large_primes, large_primes[1:]):
                temp2 *= table[q2-q1]
                Q *= temp2-1
                Q %= n
            g = gcd(Q, n)
            if g not in [1, n]:
                break
        else:
            return [(g, 1), (n//g, 1)]
    else:
        return [(n, 1)]

    return [(g, 1), (n//g, 1)]


