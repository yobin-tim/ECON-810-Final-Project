# Load packages
using Parameters, Statistics, Distributions, ProgressBars, SharedArrays, Distributed, StatsPlots
using CSV, DataFrames

# Load data
data = CSV.read("data/data_filt.csv", DataFrame)

# @df combine( grou )

data_sub = combine( groupby(data, [:age, :education_resp_HS]) , :income => median)


@df data_sub plot(:age, :income_median, group = :education_resp_HS, linewidth = 2.5, legend = :topleft, 
                xlabel = "Age", ylabel = "Median Income ( '000s dollars )", 
                title = "Median Income by Age and Education Level")
savefig("figures/fig_median_income_by_age_and_education.pdf")

addprocs(4) 



@everywhere include("functions.jl")
include("aux_functions.jl")

path = "Saved Results/"

begin 
    # @time prim, res, pre_comp = Init();
    # @time vfn(prim, res, pre_comp)
    # save_data(res, path, "_v2")

    # prim, res, pre_comp = Init()
    # vfn(prim,res,pre_comp)

    # prim, res_vmu, pre_comp = Init(μarg = 0.95)
    # vfn(prim,res_vmu,pre_comp)
    # save_data(res, path, "_vmu")

    # prim, res_vA, pre_comp = Init(Aarg = 1.05)
    # vfn(prim,res_vA,pre_comp)
    # save_data(res, path, "_vA")


    # prim, res_vthr, pre_comp = Init(thr = 0.5)
    # vfn(prim,res_vthr,pre_comp)
    # save_data(res, path)
end 

prim, res, pre_comp = read_data(path, "_valt")

sim = Init2(prim)
runsim(prim, res, sim, pre_comp, 4.0)

data_sim= sim_to_df(sim, prim)
data_sim = data_sim[data_sim.year .< 30, :]
data_sim[!, :college] = data_sim[!, :S] .> 4.0
data_sim[!, :emp_simple] = data_sim.emp_status .== 0.0

@df combine( groupby(data_sim, :year) , :emp_simple => mean) plot(:year, :emp_simple_mean)

function label(x)
    if x == 0 
        return "Unemployed"
    elseif x == 1
        return "Employed Low Firm"
    else
        return "Employed High Firm"
    end
end

data_sim[!, :emp_status_label] = label.(data_sim[:, :emp_status])


data_sim_sub = combine( groupby(data_sim, [:year, :emp_status_label]) , :hc => mean)

@df data_sim_sub plot(:year, :hc_mean, group = :emp_status_label, linewidth = 2.5, legend = :topleft, 
                xlabel = "Age", ylabel = "Median Income ", 
                title = "Median Income by (model) Age and Employment")

# savefig("figures/fig_median_income_by_age_and_employment_sim_data.pdf")

data_sim_sub = combine( groupby(data_sim, [:year, :college]) , :income => median)

@df data_sim_sub plot(:year, :income_median, group = :college, linewidth = 2.5, legend = :topleft, 
                xlabel = "Age", ylabel = "Median Income ( '000s dollars )", 
                title = "Median Income by Age and Education Level")



begin
    @df combine(groupby(data_sim[data_sim.emp_status .== 0, :], :year), :S => mean) plot(:year, :S_mean, linewidth = 2.5, legend = :topleft, 
                    xlabel = "Age", ylabel = "Median Schooling", label = "Unemployed",
                    title = "Median Schooling by Age and Employement Status")

    @df combine(groupby(data_sim[data_sim.emp_status .== 1, :], :year), :S => mean) plot!(:year, :S_mean, linewidth = 2.5, label = "Employed Low Firm")
    @df combine(groupby(data_sim[data_sim.emp_status .== 2, :], :year), :S => mean) plot!(:year, :S_mean, linewidth = 2.5, label = "Employed High Firm")                
end