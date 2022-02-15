from math import isqrt
from utils import valuation

def factor_brute(n:int) -> list[tuple[int, int]]:
    for i in range(2, isqrt(n)+1):
        if n%i == 0:
            v = valuation(n, i)
            return [(i, v)] + factor_brute(n//(i**v))
    return [(n, 1)]
