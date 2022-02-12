from math import isqrt
from utils import valuation

# Express n as difference of squares
# Used mainly for semiprimes with factors close to each other
# Complexity = 0(|p-q|) where n = p*q
def factor(n:int) ->list[tuple]:
    facs = []
    if n%2 == 0:
        k = valuation(n, 2)
        n//=(2**k)
        facs += [(2, k)]
    a = isqrt(n) + 1
    b2 = a*a - n
    while isqrt(b2)**2 != b2:
        b2+= 2*a + 1
        a+=1
    return facs + [(a-isqrt(b2),1), (a+isqrt(b2),1)]
