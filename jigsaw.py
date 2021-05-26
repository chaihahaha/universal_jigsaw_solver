from ortools.linear_solver import pywraplp
from utils import *
from data import *
import math as m


solver = pywraplp.Solver.CreateSolver('SCIP')
h = len(I0)
w = len(I0[0])
di_step = 1
dj_step = 1
dθ_step = m.pi*0.25

T = {}
for di in range(-h//2, h//2, di_step):
    for dj in range(-w//2, w//2, dj_step):
        for dθ in [n*dθ_step for n in range(int(m.pi*2/dθ_step))]:
            T[di, dj, dθ] = affine_transform_matrix(di, dj, dθ, h, w)

var = []
Ms = []
for k in range(len(Is)):
    I = Is[k]
    xTs = []
    xs = []
    for di in range(-h//2, h//2, di_step):
        for dj in range(-w//2, w//2, dj_step):
            for dθ in [n*dθ_step for n in range(int(m.pi*2/dθ_step))]:
                x = solver.IntVar(0.0, 1.0, str(k)+','+str(di)+','+str(dj)+','+"{:.2f}".format(dθ))
                xs.append(x)
                var.append(x)
                xT = mul(x, T[di, dj, dθ])
                xTs.append(xT)
    M = column(I)
    M = matmul(matsum(xTs), M)
    Ms.append(M)

    solver.Add(solver.Sum(xs) == 1)
del T

I_col = matsum(Ms)
lesseq(I_col, 1, solver)

I_mat = reshape(I_col, h, w)
J = condition_sum(mask, I_mat)
solver.Maximize(J)
status = solver.Solve()
if status == pywraplp.Solver.OPTIMAL:
    print('Solution:')
    print('Objective value =', solver.Objective().Value())
    for x in var:
        if x.solution_value() > 0:
            I_i, di, dj, dθ = str(x).split(',')
            print("piece:",I_i, "\tΔx:", dj, "\tΔy:",di, "\tΔθ", dθ)

    
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
