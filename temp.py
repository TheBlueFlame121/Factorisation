from math import sqrt

def r(b, a):
    return b%(2*a)

def rho(a, b, c):
    temp = r(-b, c)
    D = b**2 - 4*a*c
    return c, temp, (temp**2 - D)//4c

def is_reduced(a, b, c):
    sqD = sqrt(b**2 - 4*a*c)
    if abs(sqD - 2*abs(a)) < b and b < sqD:
        return True
    return False

def reduce_quadratic_form(a, b, c):
    while not is_reduced(a, b, c):
        a, b, c = rho(a, b, c)
    return a, b, c

