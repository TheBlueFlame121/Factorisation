def valuation(n:int, i:int) -> int:
    v = 0
    while n%i == 0:
        v+=1
        n//=i
    return v
