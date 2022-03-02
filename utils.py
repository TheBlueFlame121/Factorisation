from math import isqrt
from random import randint

def valuation(n:int, i:int) -> int:
    v = 0
    while n%i == 0:
        v+=1
        n//=i
    return v

def legendre_symbol(n:int, p:int) -> int:
    if n%p == 0:
        return 0
    elif pow(n, (p-1)//2, p) == 1:
        return 1
    else:
        return -1

def is_square(n:int) -> bool:
    return isqrt(n)**2 == n

def tonelli_shanks(n:int, p:int) -> tuple[int, int]:
    assert legendre_symbol(n, p) == 1, "N is not a quadratic residue"
    pm1 = p-1
    S = valuation(pm1, 2)
    Q = pm1//(2**S)
    z = randint(1, p)
    while legendre_symbol(z, p) != -1:
        z = randint(1, p)
    M = S
    c = pow(z, Q, p)
    t = pow(n, Q, p)
    R = pow(n, (Q+1)//2, p)
    while True:
        if t==0:
            return 0, 0
        if t == 1:
            return R, p-R
        i = 1
        tsq = pow(t, 2, p)
        while tsq != 1:
            tsq = pow(tsq, 2, p)
            i+=1
        b = pow(c, pow(2, M-i-1, p-1), p)
        R= (R*b)%p
        t = (t*b*b)%p
        c = (b*b)%p
        M = i

