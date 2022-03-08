using Plots, Distributed, NPZ, SharedArrays
addprocs(4)
@everywhere include("final_project_functions.jl")
@everywhere include("tauchen.jl")
prim,res = Init()


vfn(prim,res)

begin
    npzwrite("Saved Results/U_v2.npz", res.val_fun.U)
    npzwrite("Saved Results/W_L_v2.npz", res.val_fun.W_L)
    npzwrite("Saved Results/W_H_v2.npz", res.val_fun.W_H)
    npzwrite("Saved Results/k_pol_U_v2.npz", res.k_pol.U)
    npzwrite("Saved Results/s_pol_U_v2.npz", res.s_pol.U)
    # npzwrite("Saved Results/k_pol_ind_U1.npz", res.k_pol.ind_U)
    # npzwrite("Saved Results/s_pol_ind_U1.npz", res.s_pol.ind_U)
    npzwrite("Saved Results/k_pol_WL_v2.npz", res.k_pol.W_L)
    npzwrite("Saved Results/s_pol_WL_v2.npz", res.s_pol.W_L)
    # npzwrite("Saved Results/k_pol_ind_WL1.npz", res.k_pol_ind_W_L)
    # npzwrite("Saved Results/s_pol_ind_WL1.npz", res.s_pol_ind_W_L)
    npzwrite("Saved Results/k_pol_WH_v2.npz", res.k_pol.W_H)
    npzwrite("Saved Results/s_pol_WH_v2.npz", res.s_pol.W_H)
    # npzwrite("Saved Results/k_pol_ind_WH1.npz", res.k_pol_ind_W_H)
    # npzwrite("Saved Results/s_pol_ind_WH1.npz", res.s_pol_ind_W_H)
end

