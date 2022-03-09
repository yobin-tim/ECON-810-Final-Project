# Load packages
using Parameters, Statistics, Distributions, ProgressBars, SharedArrays, Distributed
using Plots, JLD, NPZ

addprocs(6) 

# Plot style
begin
    theme(:juno); #:vibrant
    default(fontfamily="Computer Modern", framestyle=:box); # LaTex-style
    gr(size = (400, 400)); # default plot size
end
res_v2 = load("Saved Results/results_v2.jld")

@everywhere include("functions.jl")


# @time prim, res, pre_comp = Init();
# @time vfn(prim, res, pre_comp)

# using JLD

# jldopen("results.jld", "w") do file
#     write(file, "val_fun", res.val_fun)
#     write(file, "k_pol", res.k_pol)
#     write(file, "s_pol", res.s_pol)
# end
prim, res, pre_comp = Init()
vfn(prim,res,pre_comp)
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

prim, res_vmu, pre_comp = Init(μarg = 0.95)
vfn(prim,res_vmu,pre_comp)
begin
    npzwrite("Saved Results/U_vmu.npz", res_vmu.val_fun.U)
    npzwrite("Saved Results/W_L_vmu.npz", res_vmu.val_fun.W_L)
    npzwrite("Saved Results/W_H_vmu.npz", res_vmu.val_fun.W_H)
    npzwrite("Saved Results/k_pol_U_vmu.npz", res_vmu.k_pol.U)
    npzwrite("Saved Results/s_pol_U_vmu.npz", res_vmu.s_pol.U)
    # npzwrite("Saved Results/k_pol_ind_vmu.npz", res_vmu.k_pol.ind_U)
    # npzwrite("Saved Results/s_pol_ind_vmu.npz", res_vmu.s_pol.ind_U)
    npzwrite("Saved Results/k_pol_WL_vmu.npz", res_vmu.k_pol.W_L)
    npzwrite("Saved Results/s_pol_WL_vmu.npz", res_vmu.s_pol.W_L)
    # npzwrite("Saved Results/k_pol_ind_vmu.npz", res_vmu.k_pol_ind_W_L)
    # npzwrite("Saved Results/s_pol_ind_vmu.npz", res_vmu.s_pol_ind_W_L)
    npzwrite("Saved Results/k_pol_WH_vmu.npz", res_vmu.k_pol.W_H)
    npzwrite("Saved Results/s_pol_WH_vmu.npz", res_vmu.s_pol.W_H)
    # npzwrite("Saved Results/k_pol_ind_vmu.npz", res_vmu.k_pol_ind_W_H)
    # npzwrite("Saved Results/s_pol_ind_vmu.npz", res_vmu.s_pol_ind_W_H)
end

prim, res_vA, pre_comp = Init(Aarg = 1.05)
vfn(prim,res_vA,pre_comp)
begin
    npzwrite("Saved Results/U_vA.npz", res_vA.val_fun.U)
    npzwrite("Saved Results/W_L_vA.npz", res_vA.val_fun.W_L)
    npzwrite("Saved Results/W_H_vA.npz", res_vA.val_fun.W_H)
    npzwrite("Saved Results/k_pol_U_vA.npz", res_vA.k_pol.U)
    npzwrite("Saved Results/s_pol_U_vA.npz", res_vA.s_pol.U)
    # npzwrite("Saved Results/k_pol_ind_vA.npz", res_vA.k_pol.ind_U)
    # npzwrite("Saved Results/s_pol_ind_vA.npz", res_vA.s_pol.ind_U)
    npzwrite("Saved Results/k_pol_WL_vA.npz", res_vA.k_pol.W_L)
    npzwrite("Saved Results/s_pol_WL_vA.npz", res_vA.s_pol.W_L)
    # npzwrite("Saved Results/k_pol_ind_vA.npz", res_vA.k_pol_ind_W_L)
    # npzwrite("Saved Results/s_pol_ind_vA.npz", res_vA.s_pol_ind_W_L)
    npzwrite("Saved Results/k_pol_WH_vA.npz", res_vA.k_pol.W_H)
    npzwrite("Saved Results/s_pol_WH_vA.npz", res_vA.s_pol.W_H)
    # npzwrite("Saved Results/k_pol_ind_vA.npz", res_vA.k_pol_ind_W_H)
    # npzwrite("Saved Results/s_pol_ind_vA.npz", res_vA.s_pol_ind_W_H)
end


prim, res_vthr, pre_comp = Init(thr = 0.5)
vfn(prim,res_vthr,pre_comp)
begin
    npzwrite("Saved Results/U_vthr.npz", res_vthr.val_fun.U)
    npzwrite("Saved Results/W_L_vthr.npz", res_vthr.val_fun.W_L)
    npzwrite("Saved Results/W_H_vthr.npz", res_vthr.val_fun.W_H)
    npzwrite("Saved Results/k_pol_U_vthr.npz", res_vthr.k_pol.U)
    npzwrite("Saved Results/s_pol_U_vthr.npz", res_vthr.s_pol.U)
    # npzwrite("Saved Results/k_pol_ind_vthr.npz", res_vthr.k_pol.ind_U)
    # npzwrite("Saved Results/s_pol_ind_vthr.npz", res_vthr.s_pol.ind_U)
    npzwrite("Saved Results/k_pol_WL_vthr.npz", res_vthr.k_pol.W_L)
    npzwrite("Saved Results/s_pol_WL_vthr.npz", res_vthr.s_pol.W_L)
    # npzwrite("Saved Results/k_pol_ind_vthr.npz", res_vthr.k_pol_ind_W_L)
    # npzwrite("Saved Results/s_pol_ind_vthr.npz", res_vthr.s_pol_ind_W_L)
    npzwrite("Saved Results/k_pol_WH_vthr.npz", res_vthr.k_pol.W_H)
    npzwrite("Saved Results/s_pol_WH_vthr.npz", res_vthr.s_pol.W_H)
    # npzwrite("Saved Results/k_pol_ind_vthr.npz", res_vthr.k_pol_ind_W_H)
    # npzwrite("Saved Results/s_pol_ind_vthr.npz", res_vthr.s_pol_ind_W_H)
end

begin
    res.val_fun.U = npzread("Saved Results/U_v2.npz")
    res.val_fun.W_L = npzread("Saved Results/W_L_v2.npz")
    res.val_fun.W_H = npzread("Saved Results/W_H_v2.npz")
    res.k_pol.U = npzread("Saved Results/k_pol_U_v2.npz"  )
    res.s_pol.U = npzread("Saved Results/s_pol_U_v2.npz"  )
    # res.k_pol.ind_U = npzread("Saved Results/k_pol_ind_U1.npz")
    # res.s_pol.ind_U = npzread("Saved Results/s_pol_ind_U1.npz")
    res.k_pol.W_L = npzread("Saved Results/k_pol_WL_v2.npz")
    res.s_pol.W_L = npzread("Saved Results/s_pol_WL_v2.npz")
    # res.k_pol_ind_W_L = npzread("Saved Results/k_pol_ind_WL1.npz")
    # res.s_pol_ind_W_L = npzread("Saved Results/s_pol_ind_WL1.npz")
    res.k_pol.W_H = npzread("Saved Results/k_pol_WH_v2.npz")
    res.s_pol.W_H = npzread("Saved Results/s_pol_WH_v2.npz")
    # res.k_pol_ind_W_H = npzread("Saved Results/k_pol_ind_WH1.npz")
    # res.s_pol_ind_W_H = npzread("Saved Results/s_pol_ind_WH1.npz")
end

sim = Init2(prim)
runsim(prim, res, sim, pre_comp)
# res_v2 = Results(res_v2["val_fun"], res_v2["k_pol"], res_v2["s_pol"])
res_v2 = res

plot(res_v2.val_fun.U[:, 50, 56, 12])
plot!(res_v2["val_fun"].W_H[:, 50, 56, 12])
plot!(res_v2["val_fun"].W_L[:, 50, 56, 12])

plot(res_v2["s_pol"].U[20,40,1,:], label="1")
plot!(res_v2["s_pol"].U[20,40,10,:], label="10")
plot!(res_v2["s_pol"].U[20,40,20,:], label="20")
plot!(res_v2["s_pol"].U[20,40,30,:], label="30")
plot(res_v2["s_pol"].U[20,40,100,:], label="100")

res_v2["s_pol"].U[:,:,1,:]