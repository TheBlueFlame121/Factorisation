import numpy as np
from numba import njit

def fast_gauss(matrix):
    matrix%=2
    m, n = matrix.shape
    mark = [False]*m

    for j, column in enumerate(matrix.T):
        #  try:
           #  pivot = np.where(column == 1)[0][0]
        i=0
        while i < n  and matrix[i, j] != 1:
           i+=1
        if i==n:
            continue
        pivot = i
        mark[pivot] = True

        for k, col in enumerate(matrix.T):
            if k == j:
                continue
            if col[pivot] == 1:
                matrix[:, k] += matrix[:, j]
                matrix[:, k] %= 2
        #  except (ValueError, IndexError):
            #  pass

    return mark, matrix

def find_perfect_square(matrix):
    print("\n[+] Running gauss elimination")
    mark, m = fast_gauss(matrix)
    #  print(m)
    m = m.tolist()
    

    independent = [(r.index(1), i,r) for i, r in enumerate(m) if mark[i]]
    dependent = [(i,r) for i, r in enumerate(m) if not mark[i]]

    for index, row in dependent:
        res = [index]

        check = 0
        if all(i==0 for i in row):
            continue

        for index, element in enumerate(row):
            if element == 1:
                for one, i, r in independent:
                    if one == index:
                        res.append(i)
                        break
                else:
                    check = 1
                    break
        if check == 1:
            continue
        yield res            
