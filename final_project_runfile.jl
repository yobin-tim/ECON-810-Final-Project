using Plots, Distributed, NPZ, SharedArrays
addprocs(4)
@everywhere include("final_project_functions.jl")
@everywhere include("tauchen.jl")
prim,res = Init()


vfn(prim,res)