using Plots, Distributed, NPZ, SharedArrays
addprocs(4)
@everywhere include("PS4_functions.jl")
@everywhere include("../Problem Set 2/tauchen.jl")
prim,res = Init()


vfn2(prim,res)
npzwrite("Problem Set 4/Saved Results/V.npz", res.V)
npzwrite("Problem Set 4/Saved Results/k_pol.npz", res.k_pol)
npzwrite("Problem Set 4/Saved Results/s_pol.npz", res.s_pol)
npzwrite("Problem Set 4/Saved Results/k_pol_ind.npz", res.k_pol_ind)
npzwrite("Problem Set 4/Saved Results/s_pol_ind.npz", res.s_pol_ind)

res.V = npzread("Problem Set 4/Saved Results/V.npz")
res.k_pol = npzread("Problem Set 4/Saved Results/k_pol.npz")
res.s_pol = npzread("Problem Set 4/Saved Results/s_pol.npz")
res.k_pol_ind = npzread("Problem Set 4/Saved Results/k_pol_ind.npz")
res.s_pol_ind = npzread("Problem Set 4/Saved Results/s_pol_ind.npz")

sim = Init2(prim)
runsim(prim, res, sim)