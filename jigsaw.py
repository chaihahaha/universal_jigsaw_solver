from ortools.linear_solver import pywraplp
from utils import *
from data import *
import math as m


solver = pywraplp.Solver.CreateSolver('SCIP')
h = len(I0)
w = len(I0[0])

var = []
Ms = []
for k in range(len(Is)):
    I = Is[k]
    xTs = []
    xs = []
    for di in range(-h//2,h//2):
        for dj in range(-w//2,w//2):
            for dθ in [0, m.pi*0.5, m.pi, m.pi*1.5]:
                T = affine_transform_matrix(di, dj, dθ, h, w)
                x = solver.IntVar(0.0, 1.0, str(k)+','+str(di)+','+str(dj)+','+str(int(dθ/(m.pi*0.5))))
                xs.append(x)
                var.append(x)
                xT = mul(x, T)
                xTs.append(xT)
    M = column(I)
    M = matmul(matsum(xTs), M)
    Ms.append(M)

    solver.Add(solver.Sum(xs) == 1)

I_col = matsum(Ms)
lessthan(I_col, 1, solver)

I_mat = reshape(I_col, h, w)
J = condition_sum(mask, I_mat)
solver.Maximize(J)
status = solver.Solve()
if status == pywraplp.Solver.OPTIMAL:
    print('Solution:')
    print('Objective value =', solver.Objective().Value())
    for x in var:
        if x.solution_value() > 0:
            I_i, di, dj, dθ = tuple(map(int, str(x).split(',')))
            print("piece:",I_i, "Δx:", dj, "Δy:",di, "Δθ", dθ)

    
    for k in range(len(Is)):
        I_ans = reshape(Ms[k],h,w)
        for i in range(h):
            for j in range(w):
                print(int(I_ans[i][j].solution_value()),end=' ')
            print()
        print()
    print("Recovered image:")
    for i in range(h):
        for j in range(w):
            print(int(I_mat[i][j].solution_value()),end=' ')
        print()
