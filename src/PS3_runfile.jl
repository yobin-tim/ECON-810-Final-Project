using StatFiles, DataFrames, CSV, Plots, Distributed, NPZ, SharedArrays, GLM

include("PS3_functions.jl")

prim,res = Init()

vfn(prim,res)
npzwrite("Problem Set 3/Saved Results/U_asset_pol.npz", res.U_asset_pol)
npzwrite("Problem Set 3/Saved Results/W_asset_pol.npz", res.W_asset_pol)
npzwrite("Problem Set 3/Saved Results/theta.npz", res.θ)
npzwrite("Problem Set 3/Saved Results/wage_rate.npz", res.wage_rate)

sim = Init2(prim)

runsim(prim,res,sim)


# theme(:vibrant) # :dark, :light, :plain, :grid, :tufte, :presentation, :none
# default(fontfamily="Computer Modern", framestyle=:box) # LaTex-style
# gr(size = (800, 600))
#Question 3
#Part A: distribution of assets from simulated data
asset_last_period   = histogram(sim.b[:,119], normalize = :probability, ylabel = "% of population", xlabel = "Asset level (start of last period)", bins = :sturges)
assets_all_periods  = histogram(reshape(sim.b,:,1), normalize = :probability, ylabel = "% of population", xlabel = "Asset level (all periods)", bins = :sturges)
plot(asset_last_period, assets_all_periods, legend = false, title = "Distribution of assets")
savefig("Problem Set 3/Figures/distribution of assets.pdf")
#Part B: distribution of wages from simulated data
wage_last_period   = histogram(sim.search_wage[:,119], normalize = :probability, ylabel = "% of population", xlabel = "Wages (start of last period)", bins = :sturges)
wage_period_40   = histogram(sim.search_wage[:,40], normalize = :probability, ylabel = "% of population", xlabel = "Wages (start of period 40)", bins = :sturges)
wage_period_80   = histogram(sim.search_wage[:,80], normalize = :probability, ylabel = "% of population", xlabel = "Wages (start of period 80)", bins = :sturges)
wage_all_periods  = histogram(reshape(sim.search_wage,:,1), normalize = :probability, ylabel = "% of population", xlabel = "Wage (all periods)", bins = :sturges)
plot(wage_last_period, wage_period_40, wage_period_80, wage_all_periods, layout  = (2,2), legend = false)
savefig("Problem Set 3/Figures/distribution of wages.pdf")
#Part C: unemployment rate
plot(1:prim.T, 1 .- mean(sim.emp_status, dims=1)', 
    title = "Unemployment rate by period",
    xlabel = "Periods", 
    xticks = 0:10:prim.T,
    legend = false)
savefig("Problem Set 3/Figures/unemployment rate.pdf")
# Part D: plot average earnings and assets over the lifecycle
plot(1:prim.T, mean(sim.income, dims=1)', label = "Earnings")
plot!(1:prim.T, mean(sim.b, dims=1)', label = "Assets")
savefig("Problem Set 3/Figures/earnings and assets.pdf")
# Part E: average gain in earnings when employed
pct_gain_earnings = wage_growth(sim)
avg_gain_earnings = mean(pct_gain_earnings)
#Part F: earnings around job loss
earnings_path = earnings_loss(sim)
avg_earnings_path = mean(earnings_path, dims=1)
earnings_path_1 = plot( -4:8, avg_earnings_path'[1:end-1],
                xlims  = (-4,8),
                xticks = -4:1:8,
                label  = "Avg. earnings",
                xlabel = "Months from layoff"
                )
savefig("Problem Set 3/Figures/avg_earnings_path.pdf")
#Part G: consumption around job loss
c_path = consumption_loss(sim)
avg_c_path = mean(c_path, dims=1)
c_path_1 = plot( -4:8, avg_c_path'[1:end-1],
                xlims  = (-4,8),
                xticks = -4:1:8,
                label  = "Avg. consumption",
                xlabel = "Months from layoff"
                )
savefig("Problem Set 3/Figures/avg_consumption_path.pdf")

###### Increase in transfer to the unemployed
prim2,res2 = Init(zarg = 0.22)
vfn(prim2,res2)
sim2 = Init2(prim2)
runsim(prim2,res2,sim2)

pct_gain_earnings_2 = wage_growth(sim2)
avg_gain_earnings_2 = mean(pct_gain_earnings_2)

earnings_path_2 = earnings_loss(sim2)
avg_earnings_path_2 = mean(earnings_path_2, dims=1)
earnings_path_2 = plot( -4:8, avg_earnings_path_2'[1:end-1],
                xlims  = (-4,8),
                xticks = -4:1:8,
                label  = "Avg. earnings (10% increase in z)",
                xlabel = "Months from layoff"
                )
c_path_2 = consumption_loss(sim2)
avg_c_path_2 = mean(c_path_2, dims=1)
c_path_2 = plot( -4:8, avg_c_path_2'[1:end-1],
                xlims  = (-4,8),
                xticks = -4:1:8,
                label  = "Avg. consumption (10% increase in z)",
                xlabel = "Months from layoff"
                )
plot(earnings_path_1,earnings_path_2, c_path_1, c_path_2, layout = (2,2))
plot(-4:8, [avg_earnings_path'[1:end-1] avg_earnings_path_2'[1:end-1]],
                xlims  = (-4,8),
                xticks = -4:1:8,
                label  = ["Avg. earnings" "Avg. earnings (10% increase in z)"],
                xlabel = "Months from layoff"
                )
savefig("Problem Set 3/Figures/avg_earnings_path_comparison.pdf")

plot(-4:8, [avg_c_path'[1:end-1] avg_c_path_2'[1:end-1]],
                xlims  = (-4,8),
                xticks = -4:1:8,
                label  = ["Avg. consumption" "Avg. consumption (10% increase in z)"],
                xlabel = "Months from layoff"
                )
savefig("Problem Set 3/Figures/avg_consumption_path_comparison.pdf")

