# Load packages
using Parameters, Statistics, Distributions, ProgressBars, SharedArrays, Distributed
using Plots, JLD

# Plot style
begin
    theme(:juno); #:vibrant
    default(fontfamily="Computer Modern", framestyle=:box); # LaTex-style
    gr(size = (400, 400)); # default plot size
end

include("functions.jl")

# addprocs(6) 

# @time prim, res, pre_comp = Init();
# @time vfn(prim, res, pre_comp)

# using JLD

# jldopen("results.jld", "w") do file
#     write(file, "val_fun", res.val_fun)
#     write(file, "k_pol", res.k_pol)
#     write(file, "s_pol", res.s_pol)
# end

res_v2 = load("Saved Results/results_v2.jld")

prim, res, pre_comp = Init()

# res_v2 = Results(res_v2["val_fun"], res_v2["k_pol"], res_v2["s_pol"])

plot(res_v2["val_fun"].U[:, 50, 56, 12])
plot!(res_v2["val_fun"].W_H[:, 50, 56, 12])
plot!(res_v2["val_fun"].W_L[:, 50, 56, 12])

plot(res_v2["s_pol"].U[20,40,1,:], label="1")
plot!(res_v2["s_pol"].U[20,40,10,:], label="10")
plot!(res_v2["s_pol"].U[20,40,20,:], label="20")
plot!(res_v2["s_pol"].U[20,40,30,:], label="30")
plot(res_v2["s_pol"].U[20,40,100,:], label="100")

res_v2["s_pol"].U[:,:,1,:]