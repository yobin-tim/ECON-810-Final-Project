# Load packages
using Parameters, Statistics, Distributions, ProgressBars, SharedArrays, Distributed
addprocs(4) 

@everywhere include("functions.jl")
include("aux_functions.jl")

path = "Saved Results/"



begin 
    # @time prim, res, pre_comp = Init();
    # @time vfn(prim, res, pre_comp)
    # save_data(res, path, "_v2")

    # prim, res_vthr, pre_comp = Init(thr = (31/151))
    # vfn(prim,res_vthr,pre_comp)
    # save_data(res_vthr, path, "_vthr")

    prim, res, pre_comp = Init()
    vfn(prim,res,pre_comp)
    save_data(res, path, "_vbig")

    # prim, res_vmu, pre_comp = Init(μarg = 0.95)
    # vfn(prim,res_vmu,pre_comp)
    # save_data(res_vmu, path, "_vmu")

    # prim, res_vA, pre_comp = Init(Aarg = 1.05)
    # vfn(prim,res_vA,pre_comp)
    # save_data(res_vA, path, "_vA")


    
end 

#Policy function graphs
pf1 =  plot(1:prim.T, [prim.s_grid[res.s_pol.U[21,50,5,:]]  prim.s_grid[res.s_pol.U[21,50,21,:]]], xlabel = "Model Age", label = ["s (low S)" "s (high S)"], title = "s (for low Sbar)")
pf2 =  plot(1:prim.T, [prim.s_grid[res_vthr.s_pol.U[21,50,5,:]]  prim.s_grid[res_vthr.s_pol.U[21,50,21,:]]], xlabel = "Model Age", label = ["s (low S)" "s (high S)"], title = "s (for high Sbar)")
plot(pf1, pf2, layout = (1,2))
savefig("figures/s_pol_U.pdf")

prim, res, pre_comp = read_data(path, "_v1")

sim = Init2(prim)
runsim(prim, res, sim, pre_comp)

# plot(res.val_fun.U[:, 50, 56, 12])
# plot!(res.val_fun.W_H[:, 50, 56, 12])
# plot!(res.val_fun.W_L[:, 50, 56, 12])

# plot(res_v2["s_pol"].U[20,40,1,:], label="1")
# plot!(res_v2["s_pol"].U[20,40,10,:], label="10")
# plot!(res_v2["s_pol"].U[20,40,20,:], label="20")
# plot!(res_v2["s_pol"].U[20,40,30,:], label="30")
# plot(res_v2["s_pol"].U[20,40,100,:], label="100")

plot(mean(prim.S_grid[sim.S], dims = 1)', title = "Mean years of Education", xlabel = "Model Time", legend = false)
plot(mean(prim.s_grid[sim.s], dims = 1)', title = "Investment in education each period", xlabel = "Model Time", legend = false)

# res_v2["s_pol"].U[:,:,1,:]
plot( mean( prim.h_grid[sim.hc] , dims = 1)' , legend = false)
plot( mean( prim.S_grid[sim.S] , dims = 1)' , legend = false)
plot( mean( prim.s_grid[sim.s] , dims = 1)' , legend = false)

plot(100 * sum(prim.S_grid[sim.S]  .>= prim.S_bar , dims = 1)' / prim.nSim)

using StatsPlots

histogram(prim.k_grid[sim.k][:,1], bins = 50)
histogram!(prim.k_grid[sim.k][:,end-10], bins = 50, alpha=0.5)