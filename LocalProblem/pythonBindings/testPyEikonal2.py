import numpy as np
import pyEikonal2 as eik
tria=[[0.2,0.3],
       [0.,1.],
       [1.,0.],
       ]
values=np.array([[0.],[0]],dtype=np.float64)
M=np.array([
    [9.0,0.0],
    [0.0,9.0]]
    )
#sol=eik.solveLocalProblemIsotropic(tetra,values)
oopt=eik.OptimizationOptions(); # if you want to change defualts
lsopt=eik.LineSearchOptions(); # if you want to change defaults
lsopt.maxIter=10; #max iterations for backtracking
oopt.maxIter=12; #max iterations for the solver
eik.setLineSearchOptions(lsopt)
eik.setOptimizationOptions(oopt)
sol=eik.solveLocalProblem(tria,values,M)
print(sol.value)
print(sol.l)
print(sol.status)


