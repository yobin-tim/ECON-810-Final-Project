using Plots, Distributed, NPZ, SharedArrays
addprocs(4)
@everywhere include("final_project_functions.jl")
@everywhere include("tauchen.jl")
prim,res = Init()


vfn(prim,res)

begin
    npzwrite("Saved Results/U1.npz", res.U)
    npzwrite("Saved Results/W_L1.npz", res.W_L)
    npzwrite("Saved Results/W_H1.npz", res.W_H)
    npzwrite("Saved Results/k_pol_U1.npz", res.k_pol_U)
    npzwrite("Saved Results/s_pol_U1.npz", res.s_pol_U)
    npzwrite("Saved Results/k_pol_ind_U1.npz", res.k_pol_ind_U)
    npzwrite("Saved Results/s_pol_ind_U1.npz", res.s_pol_ind_U)
    npzwrite("Saved Results/k_pol_WL1.npz", res.k_pol_W_L)
    npzwrite("Saved Results/s_pol_WL1.npz", res.s_pol_W_L)
    npzwrite("Saved Results/k_pol_ind_WL1.npz", res.k_pol_ind_W_L)
    npzwrite("Saved Results/s_pol_ind_WL1.npz", res.s_pol_ind_W_L)
    npzwrite("Saved Results/k_pol_WH1.npz", res.k_pol_W_H)
    npzwrite("Saved Results/s_pol_WH1.npz", res.s_pol_W_H)
    npzwrite("Saved Results/k_pol_ind_WH1.npz", res.k_pol_ind_W_H)
    npzwrite("Saved Results/s_pol_ind_WH1.npz", res.s_pol_ind_W_H)
end

