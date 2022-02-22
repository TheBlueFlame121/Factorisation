from collections import namedtuple
from random import randint
from prime_sieves import sieve_of_Eratosthenes
from math import log, gcd

Curve = namedtuple("Curve", ["A", "B", "C", "mod"])
Point = namedtuple("Point", ["X", "Z"])

inf = Point(0, 0)

def addh(P1:Point, P2:Point, Pm:Point, C:Curve) -> Point:
    if P1 == inf:
        #  return P2
        return Point(P2[0], P2[1])
    if P2 == inf:
        #  return P1
        return Point(P1[0], P1[1])

    Z1Z2 = (P1[1]*P2[1]) % C[3]

    X1X2_AZ1Z2_sq = pow((P1[0]*P2[0] - C[0]*Z1Z2), 2, C[3])
    X1Z2pX2Z1pCZ1Z2 = (P1[0]*P2[1] + P2[0]*P1[1] + C[2]*Z1Z2) % C[3]

    Xp = (Pm[1]*(X1X2_AZ1Z2_sq - 4*C[1]*X1Z2pX2Z1pCZ1Z2*Z1Z2)) % C[3]
    Zp = (Pm[0]*pow(P1[0]*P2[1] - P2[0]*P1[1], 2)) % C[3]

    return Point(Xp, Zp)


def doubleh(P1:Point, C:Curve) -> Point:
    (X1, Z1), (A, B, C, mod) = P1, C

    temp = pow((X1*X1 - A*Z1*Z1), 2, mod)
    Xp = (temp - 4*B*(2*X1 + C*Z1)*Z1*Z1*Z1) % mod
    Zp = (4*Z1*(X1*X1*(X1 + C*Z1) + Z1*Z1*(A*X1 + B*Z1))) % mod

    return Point(Xp, Zp)


def multiply(n:int, P:Point, C:Curve) -> Point:
    if n == 0:
        return Point(0, 0)
    if n == 1:
        return Point(P[0], P[1])
    if n == 2:
        return doubleh(P, C)
    UV = Point(P[0], P[1])
    TW = doubleh(P, C)
    B = n.bit_length()
    bit_string = bin(n)[2:-1][::-1]
    for b in bit_string:
        if b == '1':
            UV = addh(TW, UV, P, C)
            TW = doubleh(TW, C)
        else:
            TW = addh(UV, TW, P, C)
            UV = doubleh(UV, C)

    if bit_string[-1] == '1':
        return addh(UV, TW, P, C)
    return doubleh(UV, C)

def factor_ecm(N:int, B1:int=10000, B2:int=None, D:int=100) -> list[tuple[int, int]]:
    if not B2:
        B2 = 100*B1
    assert B1 % 2==0 and B2 % 2 == 0

    sigma = randint(6, N-1)
    u = (sigma*sigma - 5) % N
    v = (4 * sigma) % N
    C = ((v-u)**3 * (3*u+v))*pow(4*u**3*v, -1, N) - 2
    C %= N
    curve = Curve(1, 0, C, N)
    Q = Point(pow(u, 3, N), pow(v, 3, N))
    
    # Stage 1
    for p in sieve_of_Eratosthenes(B1):
        a = int(log(B1, p))
        Q = multiply(p**a, Q, curve)
    g = gcd(int(Q[1]), N)
    if g not in [1, N]:
        return [(g, 1), (N//g, 1)]

    # Stage 2
    S = [None, doubleh(Q, curve)]
    Beta = [None, (S[1][0]*S[1][1])%N]
    S.append(doubleh(S[1], curve))
    Beta.append((S[2][0]*S[2][1])%N)
    for d in range(3, D+1):
        S.append(addh(S[d-1], S[1], S[d-2], curve))
        Beta.append((S[d][0]*S[d][1]) % N)

    g = 1
    B = B1 - 1
    T = multiply(B - 2*D, Q, curve)
    R = multiply(B, Q, curve)
    primes = sieve_of_Eratosthenes(B2 + 2*D)
    for r in range(B, B2, 2*D):
        alpha = (R[0]*R[1]) % N
        for q in [p for p in primes if (p > r+1 and p< r+ 2*D + 1)]:
            delta = (q-r)//2
            g = (g*((R[0] - S[delta][0])*(R[1] + S[delta][1]) - alpha - Beta[delta])) % N
        R, T = addh(R, S[D], T, curve), R
    print("POG")
    g = gcd(g, N)
    if g not in [1, N]:
        return [(g, 1), (N//g, 1)]
