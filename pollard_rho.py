from random import randint
from math import gcd

# Pollard's rho factorisation using Brent's cycle detection method
# Assumes n is semiprime
def factor(n:int, funcs:int, iters:int):
    for _ in range(funcs):
        c = randint(1, n-1)
        F = lambda x: (x*x + c) % n
        
        power = lam = 1
        x = randint(0, n-1)
        y = F(x)
        for i in range(iters):
            if power == lam:
                x = y
                power *= 2
                lam = 0
            y = F(y)
            lam+=1
            g = gcd(abs(x-y), n)
            if g != 1:
                return [(g, 1), (n//g, 1)]

