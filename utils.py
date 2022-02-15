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
        return 0
