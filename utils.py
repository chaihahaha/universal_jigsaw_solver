import math as m
def matmul(X, Y):
    hx = len(X)
    wx = len(X[0])
    hy = len(Y)
    wy = len(Y[0])
    assert wx==hy
    return [[sum([X[i][k] * Y[k][j] for k in range(wx)]) for j in range(wy)] for i in range(hx)]

def lesseq(M, n, solver):
    h = len(M)
    w = len(M[0])
    for i in range(h):
        for j in range(w):
            solver.Add(M[i][j] <= n)
    return

def mul(k, X):
    h = len(X)
    w = len(X[0])
    return [[k*X[i][j]  for j in range(w)] for i in range(h)]

def matsum(Ms):
    M0 = Ms[0]
    h = len(M0)
    w = len(M0[0])
    for i in range(len(Ms)):
        assert h==len(Ms[i]) and w==len(Ms[i][0])
    return [[sum([M[i][j] for M in Ms]) for j in range(w)] for i in range(h)]

def column(X):
    h = len(X)
    w = len(X[0])
    return [[X[i][j]] for i in range(h) for j in range(w)]

def reshape(M, h, w):
    r = []
    for i in range(h):
        r.append([M[i*w+j][0] for j in range(w)])
    return r

def condition_sum(mask, M):
    h = len(M)
    w = len(M[0])
    x = []
    for i in range(h):
        for j in range(w):
            if mask[i][j]:
                x.append(M[i][j])

    return sum(x)

def affine_transform_matrix(di, dj, dθ, h, w):
    R = [[0]*h*w for i in range(h*w)]
    T = [[0]*h*w for i in range(h*w)]
    ci, cj = h//2, w//2
    for ii in range(h):
        for jj in range(w):
            for i in range(h):
                for j in range(w):
                    res = dθ % m.pi/2
                    tolerance = abs(m.cos(res)) if res > 0 else 1e-4
                    if abs((ii-ci) - (m.cos(dθ)*(i-ci)-m.sin(dθ)*(j-cj))) <=tolerance and abs((jj-cj) - (m.sin(dθ)*(i-ci)+m.cos(dθ)*(j-cj))) <=tolerance:
                        R[ii*w+jj][i*w+j]=1
                        break
                else:
                    continue
                break
    for i in range(h):
        for j in range(w):
            for ii in range(h):
                for jj in range(w):
                    if (ii-i)==di and (jj-j)==dj:
                        T[ii*w+jj][i*w+j]=1
    return matmul(T, R)
