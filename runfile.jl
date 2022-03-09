# Load packages
using Parameters, Statistics, Distributions, ProgressBars, SharedArrays, Distributed
addprocs(6) 

@everywhere include("functions.jl")
include("aux_functions.jl")

path = "Saved Results/"

begin 
    # @time prim, res, pre_comp = Init();
    # @time vfn(prim, res, pre_comp)
    # save_data(res, path, "_v2")
    # using JLD

    # jldopen("results.jld", "w") do file
    #     write(file, "val_fun", res.val_fun)
    #     write(file, "k_pol", res.k_pol)
    #     write(file, "s_pol", res.s_pol)
    # end
    # prim, res, pre_comp = Init()
    # vfn(prim,res,pre_comp)
    # save_data(res, path, "_v1")
    
    prim, res_vmu, pre_comp = Init(μarg = 0.95)
    vfn(prim,res_vmu,pre_comp)
    save_data(res_vmu, path, "_vmu")

        

    # prim, res_vA, pre_comp = Init(Aarg = 1.05)
    # vfn(prim,res_vA,pre_comp)
    # save_data(res, path)


    prim, res_vthr, pre_comp = Init(thr = 0.5)
    vfn(prim,res_vthr,pre_comp)
    save_data(res, path, "_vthr")
end 

# prim, res, pre_comp = read_data(path, "_vA")

save_data(res, path, "_v1")

sim = Init2(prim)
runsim(prim, res, sim, pre_comp)

# plot(res_v2.val_fun.U[:, 50, 56, 12])
# plot!(res_v2["val_fun"].W_H[:, 50, 56, 12])
# plot!(res_v2["val_fun"].W_L[:, 50, 56, 12])

# plot(res_v2["s_pol"].U[20,40,1,:], label="1")
# plot!(res_v2["s_pol"].U[20,40,10,:], label="10")
# plot!(res_v2["s_pol"].U[20,40,20,:], label="20")
# plot!(res_v2["s_pol"].U[20,40,30,:], label="30")
# plot(res_v2["s_pol"].U[20,40,100,:], label="100")

# res_v2["s_pol"].U[:,:,1,:]