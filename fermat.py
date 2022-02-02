from math import isqrt


# Express n as difference of squares
# Used mainly for semiprimes with factors close to each other
# Complexity = 0(|p-q|) where n = p*q
def basic(n:int):
    a = isqrt(n) + 1
    b2 = a*a - n
    while isqrt(b2)**2 != b2:
        b2+= 2*a + 1
        a+=1
    return [(a-isqrt(b2),1), (a+isqrt(b),1)]
