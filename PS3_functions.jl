@everywhere using Parameters, Statistics, Distributions, ProgressBars, SharedArrays, DataFrames, CategoricalArrays

@with_kw struct Primitives
    r::Float64              = 0.04
    β::Float64              = 0.99      #discount rate    
    βf::Float64             = 0.99      #discount rate for firms and lenders
    q::Float64              = 1/(1+r)
    δ::Float64              = 0.1       #job destruction rate, quarterly
    ζ::Float64              = 1.6       #matching elasticity parameter
    M::Function             = (u,v) -> (u * v)/(u^ζ + v^ζ)^(1/ζ)
    κ::Float64              = .995      #vacancy posting cost
    z::Float64              = 0.2       #transfer from the govt
    T::Int64                = 120       #number of periods
    σ::Int64                = 2         #
    u::Function             = (c) -> (c^(1 - σ)-1)/(1 - σ)
    f::Function             = (h) -> h
    p::Function             = (θ, ζ) -> (θ/(1 + θ^ζ)^(1/ζ))
    pf::Function            = (θ, ζ) -> (1/(1 + θ^ζ)^(1/ζ))
    # pf_inv::Function        = (pf,ζ) -> (pf.^(-ζ).-1).^(1/ζ)
    τ::Float64              = 0.2       #marginal tax rate
    p_hL::Float64           = 0.5       #prob of HC decline
    p_hH::Float64           = 0.05      #prob of HC increase
    #human capital grid
    h_min::Float64              = 0.5
    h_max::Float64              = 1.5
    n_hPoints::Int64            = 25
    Δ::Float64                  = (h_max - h_min)/(n_hPoints - 1)
    h_grid::Array{Float64,1}    = range(h_min, h_max, length=n_hPoints)
    #wage grid
    ω_min::Float64              = 0.2
    ω_max::Float64              = 1.0
    n_ωPoints::Int64            = 41
    ω_grid::Array{Float64,1}    = range(ω_min, ω_max, length=n_ωPoints)
    #asset grid
    b_min::Float64              = 0.0
    b_max::Float64              = 4.0
    n_bPoints::Int64            = 80
    b_grid::Array{Float64,1}    = range(b_min, b_max, length=n_bPoints)
    #simulations
    nSim                       = 10000
end

mutable struct Results
    U::SharedArray{Float64,3}
    W::SharedArray{Float64,4}
    J::SharedArray{Float64,3}
    U_asset_pol::SharedArray{Float64,3}
    W_asset_pol::SharedArray{Float64,4}
    wage_rate::SharedArray{Float64,3}
    θ::SharedArray{Float64,3}
end


function Init(;zarg=0.2)
    prim                        = Primitives(z=zarg)
    U                           = SharedArray{Float64}(prim.n_bPoints, prim.n_hPoints,prim.T)
    W                           = SharedArray{Float64}(prim.n_ωPoints,prim.n_bPoints, prim.n_hPoints,prim.T)
    J                           = SharedArray{Float64}(prim.n_ωPoints, prim.n_hPoints,prim.T)
    U_asset_pol                 = SharedArray{Float64}(prim.n_bPoints, prim.n_hPoints,prim.T)
    W_asset_pol                 = SharedArray{Float64}(prim.n_ωPoints,prim.n_bPoints, prim.n_hPoints,prim.T)
    wage_rate                   = SharedArray{Float64}(prim.n_bPoints, prim.n_hPoints,prim.T)
    θ                           = SharedArray{Float64}(prim.n_ωPoints, prim.n_hPoints,prim.T)
    res                         = Results(U,W,J,U_asset_pol,W_asset_pol,wage_rate,θ)
    return prim, res
end

function pf_inv(κ::Float64, J::Float64)
    @unpack ζ = Primitives()
    if κ > J  #free entry condition not binding implies θ = 0.
        return 0.0
    else
        ((κ/J).^(-ζ).-1).^(1/ζ)
    end
end

function vfn(prim::Primitives, res::Results)
    @unpack β, βf, δ, ζ, T, f, h_grid, b_grid, κ, q, τ, ω_grid, n_ωPoints, n_bPoints, n_hPoints, p_hL, p_hH, z, p, u = prim
    @unpack U, W, J, θ, W_asset_pol, U_asset_pol, wage_rate = res

    J[:,:,T]    = (1 .-ω_grid) * f.(h_grid)'
    θ[:,:,T]    = pf_inv.(κ,J[:,:,T])
    for ωbh in 1:(n_ωPoints * n_bPoints * n_hPoints)
        ω,b,h       = Tuple(CartesianIndices((n_ωPoints,n_bPoints,n_hPoints))[ωbh])
        W[ω,b,h,T]  = (1 - τ) * ω_grid[ω] * f.(h_grid[h]) + b_grid[b]
    end 
    U[:,:,T]    = (z .+ b_grid) .* ones(n_bPoints,n_hPoints)
    for bh in 1:(n_bPoints * n_hPoints) 
        b,h                 = Tuple(CartesianIndices((n_bPoints,n_hPoints))[bh])
        w_index             = argmax(p.(θ[:,h,120],ζ) .*  W[:,b,h,120] .+ (1 .- p.(θ[:,h,120],ζ)) .* U[b,h,120])
        wage_rate[b,h,120]    = ω_grid[w_index]
    end

    #For the vfi, note that all matrices in res have h and t common, so they will form the
    #outermost loops. 
    for t in ProgressBar(T-1:-1:1) 
        @sync @distributed for h in 1:n_hPoints #iterate over human capital
            #for the employed value function, we need to first iterate over wages
            for ω in 1:n_ωPoints #iterate over wage rates
                #firm value function
                J[ω,h,t]    =  (1 - ω_grid[ω]) * f(h_grid[h])
                J[ω,h,t]    += βf * (1 - δ) * p_hH * J[ω,min(h+1,n_hPoints),t+1]
                J[ω,h,t]    += βf * (1 - δ) * (1 - p_hH) * J[ω,h,t+1]
                #labor market tightness
                θ[ω,h,t]    = pf_inv(κ,J[ω,h,t])

                for b in 1:n_bPoints #iterate over asset holdings
                    budget_e = (1 - τ) * ω_grid[ω] * f(h_grid[h]) + b_grid[b]
                    cand_val_e = -Inf
                    for bprime in 1:n_bPoints #iterate over asset holdings tomorrow 
                        c_e = budget_e - b_grid[bprime]
                        if c_e<0 # Since all future higher values of b' are infeasible
                            break
                        end
                        val_func_e =  u(c_e)
                        val_func_e += β * p_hH * ( (1 - δ) * W[ω,bprime,min(h+1,n_hPoints),t+1] + δ * U[bprime,min(h+1,n_hPoints),t+1])
                        val_func_e += β * (1 - p_hH) * ( (1 - δ) * W[ω,bprime,h,t+1] + δ * U[bprime,h,t+1])
                        if val_func_e > cand_val_e
                            W[ω,b,h,t]              = val_func_e
                            W_asset_pol[ω,b,h,t]    = b_grid[bprime]
                            cand_val_e              = val_func_e 
                        end
                    end
                end
            end
            #for the unemployed value function, iterate over b_grid
            for b in 1:n_bPoints 
                budget_u = z + b_grid[b]
                cand_val_u = -Inf
                for bprime in 1:n_bPoints 
                    c_u = budget_u - q * b_grid[bprime]
                    if c_u<0 # Since all future higher values of b' are infeasible
                        break
                    end
                    val_func_u    =  u(c_u)
                    # exp_next_grid = p_hL .* (p.(θ[:,max(h-1,1),t+1],ζ) .* W[:,bprime,max(h-1,1),t+1] .+ (1 .- p.(θ[:,max(h-1,1),t+1],ζ)) .* U[bprime,max(h-1,1),t+1])
                    # exp_next_grid += (1 .- p_hL) .* (p.(θ[:,h,t+1],ζ) .*  W[:,bprime,h,t+1] .+ (1 .- p.(θ[:,h,t+1],ζ)) .* U[bprime,h,t+1])
                    # val_func_u    += β .* maximum(exp_next_grid)  #todo: max of expectation instead of expectation of max yields a difference?
                    v_for_hprime  = maximum(p.(θ[:,max(h-1,1),t+1],ζ) .* W[:,bprime,max(h-1,1),t+1] .+ (1 .- p.(θ[:,max(h-1,1),t+1],ζ)) .* U[bprime,max(h-1,1),t+1])
                    v_for_h       = maximum(p.(θ[:,h,t+1],ζ) .*  W[:,bprime,h,t+1] .+ (1 .- p.(θ[:,h,t+1],ζ)) .* U[bprime,h,t+1])
                    val_func_u    += p_hL .* v_for_hprime .+ (1 .- p_hL) .* v_for_h
                    if val_func_u > cand_val_u
                        U[b,h,t]              = val_func_u
                        U_asset_pol[b,h,t]    = b_grid[bprime]
                        cand_val_u            = val_func_u 
                    end
                end
            end
        end
        #solve the job search problem, find the wage rate where the agent will search
        for bh in 1:(n_bPoints * n_hPoints) 
            b,h                 = Tuple(CartesianIndices((n_bPoints,n_hPoints))[bh])
            w_index             = argmax(p.(θ[:,h,t],ζ) .*  W[:,b,h,t] .+ (1 .- p.(θ[:,h,t],ζ)) .* U[b,h,t])
            wage_rate[b,h,t]    = ω_grid[w_index]
        end
    end
end

@everywhere mutable struct Simulations
    b::Array{Float64,2}
    b_ind::Array{Int64,2}
    c::Array{Float64,2}
    income::Array{Float64,2}
    hc::Array{Float64,2}
    hc_ind::Array{Int64,2}
    emp_status::Array{Int64,2}
    emp_streak::Array{Int64,2}
    search_wage::Array{Float64,2}
    sw_ind::Array{Int64,2}
end

function Init2(prim::Primitives)
    b           = zeros(prim.nSim, prim.T)
    b_ind       = ones(prim.nSim, prim.T)
    c           = zeros(prim.nSim, prim.T)
    income      = zeros(prim.nSim, prim.T)
    hc          = zeros(prim.nSim, prim.T)
    hc_ind      = ones(prim.nSim, prim.T)
    emp_status  = zeros(prim.nSim, prim.T)
    emp_streak  = zeros(prim.nSim, prim.T)
    search_wage = zeros(prim.nSim, prim.T)
    sw_ind      = ones(prim.nSim, prim.T)
    sim         = Simulations(b, b_ind, c, income, hc, hc_ind, emp_status, emp_streak, search_wage, sw_ind)
    return sim
end

function runsim(prim::Primitives, res::Results, sim::Simulations)
    @unpack T, h_grid, nSim, p, ζ, ω_grid, δ, z, Δ, p_hL, p_hH, n_hPoints, b_grid,q, τ = prim
    @unpack U_asset_pol, W_asset_pol, θ, wage_rate = res
    @unpack b, b_ind, c, income, hc, hc_ind, emp_status, emp_streak, search_wage, sw_ind = sim

    #Initial period setup; indexes initialized at 1, assets at 0 already
    # emp_status initialized at 0 too.
    hc[:,1]     .= h_grid[1]

    for t in 1:T
        #job finding rate
        jf_rate = p.(θ[:,:,t],ζ)
        #for each agent
        if t == 1
            for i in 1:nSim
                search_wage[i,t]    = wage_rate[b_ind[i,t],hc_ind[i,t],t]
                sw_ind[i,t]         = findall(x -> x == search_wage[i,t], ω_grid)[1]
                # matches with an employer
                if rand(Uniform(0,1)) < jf_rate[sw_ind[i,t],hc_ind[i,t]]
                    income[i,t]     = (1 - τ) * search_wage[i,t] * hc[i,t]
                    b[i,t+1]        = W_asset_pol[sw_ind[i,t],b_ind[i,t],hc_ind[i,t],t]
                    emp_status[i,t] = 1
                    emp_streak[i,t] = 1
                    #human capital adjustment
                    if rand(Uniform(0,1)) < p_hH
                        hc[i,t+1]       = min(hc[i,t] + Δ,h_grid[end])
                        hc_ind[i,t+1]   = min(hc_ind[i,t]+1, n_hPoints)
                    else
                        hc[i,t+1]       = hc[i,t]
                        hc_ind[i,t+1]   = hc_ind[i,t]
                    end
                else #no match
                    income[i,t]     = z
                    b[i,t+1]        = U_asset_pol[b_ind[i,t],hc_ind[i,t],t]
                    emp_status[i,t] = 0
                    emp_streak[i,t] = 0
                    #human capital adjustment
                    if rand(Uniform(0,1)) < p_hL
                        hc[i,t+1]       = max(hc[i,t] - Δ,h_grid[1])
                        hc_ind[i,t+1]   = max(hc_ind[i,t]-1, 1)
                    else
                        hc[i,t+1]       = hc[i,t]
                        hc_ind[i,t+1]   = hc_ind[i,t]
                    end
                end
                b_ind[i,t+1]        = findall(x -> x == b[i,t+1], b_grid)[1]
                c[i,t]              = income[i,t] + b[i,t] - q * b[i,t+1]
            end
        elseif t == T #if terminal period
            for i in 1:nSim
                #Unemployment shocks realized at the beginning of the period; no shocks before period 1
                #If unemployed or at t=1
                if emp_status[i,t-1] == 0 || (emp_status[i,t-1] == 1 && rand(Uniform(0,1)) < δ)
                    #find the wage over which they search
                    search_wage[i,t] = wage_rate[b_ind[i,t],hc_ind[i,t],t]
                    sw_ind[i,t]      = findall(x -> x == search_wage[i,t], ω_grid)[1]
                    # matches with an employer
                    if rand(Uniform(0,1)) < jf_rate[sw_ind[i,t],hc_ind[i,t]]
                        income[i,t]     = (1 - τ) * search_wage[i,t] * hc[i,t]
                        emp_status[i,t] = 1
                        if emp_status[i,t-1] == 1
                            emp_streak[i,t] = emp_streak[i,t-1] + 1
                        else
                            emp_streak[i,t] = 1
                        end
                    else #no match
                        income[i,t]     = z
                        emp_status[i,t] = 0
                        emp_streak[i,t] = 0
                    end
                    c[i,t]              = income[i,t] + b[i,t]                    
                else #if employed
                    income[i,t]         = (1 - τ) * search_wage[i,t] * hc[i,t]
                    c[i,t]              = income[i,t] + b[i,t]
                    emp_status[i,t]     = 1
                    emp_streak[i,t]     = emp_streak[i,t-1] + 1
                end
            end
        else #if non-extreme periods
            for i in 1:nSim
                #Unemployment shocks realized at the beginning of the period; no shocks before period 1
                #If unemployed or at t=1
                if emp_status[i,t-1] == 0 || (emp_status[i,t-1] == 1 && rand(Uniform(0,1)) < δ)
                    #find the wage over which they search
                    search_wage[i,t]    = wage_rate[b_ind[i,t],hc_ind[i,t],t]
                    sw_ind[i,t]         = findall(x -> x == search_wage[i,t], ω_grid)[1]
                    # matches with an employer
                    if rand(Uniform(0,1)) < jf_rate[sw_ind[i,t],hc_ind[i,t]]
                        income[i,t]     = (1 - τ) * search_wage[i,t] * hc[i,t]
                        b[i,t+1]        = W_asset_pol[sw_ind[i,t],b_ind[i,t],hc_ind[i,t],t]
                        emp_status[i,t] = 1
                        if emp_status[i,t-1] == 1
                            emp_streak[i,t] = emp_streak[i,t-1] + 1
                        else
                            emp_streak[i,t] = 1
                        end
                        #human capital adjustment
                        if rand(Uniform(0,1)) < p_hH
                            hc[i,t+1]       = min(hc[i,t] + Δ,h_grid[end])
                            hc_ind[i,t+1]   = min(hc_ind[i,t]+1, n_hPoints)
                        else
                            hc[i,t+1]       = hc[i,t]
                            hc_ind[i,t+1]   = hc_ind[i,t]
                        end
                    else #no match
                        income[i,t]     = z
                        b[i,t+1]        = U_asset_pol[b_ind[i,t],hc_ind[i,t],t]
                        emp_status[i,t] = 0
                        emp_streak[i,t] = 0
                        #human capital adjustment
                        if rand(Uniform(0,1)) < p_hL
                            hc[i,t+1]       = max(hc[i,t] - Δ,h_grid[1])
                            hc_ind[i,t+1]   = max(hc_ind[i,t]-1, 1)
                        else
                            hc[i,t+1]       = hc[i,t]
                            hc_ind[i,t+1]   = hc_ind[i,t]
                        end
                    end
                    b_ind[i,t+1]        = findall(x -> x == b[i,t+1], b_grid)[1]
                    c[i,t]              = income[i,t] + b[i,t] - q * b[i,t+1]
                else #if employed
                    search_wage[i,t]    = search_wage[i,t-1]
                    sw_ind[i,t]         = sw_ind[i,t-1]
                    income[i,t]         = (1 - τ) * search_wage[i,t] * hc[i,t]
                    b[i,t+1]            = W_asset_pol[sw_ind[i,t],b_ind[i,t],hc_ind[i,t],t]
                    b_ind[i,t+1]        = findall(x -> x == b[i,t+1], b_grid)[1]
                    c[i,t]              = income[i,t] + b[i,t] - q * b[i,t+1]
                    emp_status[i,t]     = 1
                    emp_streak[i,t]     = emp_streak[i,t-1] + 1
                    #human capital adjustment
                    if rand(Uniform(0,1)) < p_hH
                        hc[i,t+1]       = min(hc[i,t] + Δ,h_grid[end])
                        hc_ind[i,t+1]   = min(hc_ind[i,t]+1, n_hPoints)
                    else
                        hc[i,t+1]       = hc[i,t]
                        hc_ind[i,t+1]   = hc_ind[i,t]
                    end
                end
            end
        end
    end
end

function wage_growth(sim::Simulations,streak::Int64=24)
    @unpack emp_streak, search_wage, income, hc = sim
    set = findall(emp_streak .>= streak)
    wage_change = zeros(size(set,1),1)
    for ind in 1:size(set,1) 
        w_final        = income[set[ind]]
        w_initial      = income[set[ind][1], set[ind][2] - streak + 1]
        wage_change[ind]    = (w_final - w_initial)/w_initial * 100
    end
    return wage_change
end

function earnings_loss(sim::Simulations,preloss::Int64=5,postloss::Int64=8)
    @unpack T = Primitives()
    @unpack emp_streak, search_wage, hc, income = sim
    matrix = emp_streak[:,preloss+1:end]  #omitting the first period where everyone is unemployed
    set = findall(iszero,matrix)
    wage_path = zeros(size(set,1),preloss+postloss+1)
    for ind in 1:size(set,1) 
        first_index         = set[ind][2]
        last_index          = min(set[ind][2] + preloss + postloss,T)
        if income[set[ind][1], set[ind][2]+ preloss-1] == 0
            continue
        end
        # if findall(isone,wage[set[ind][1], 1:set[ind][2]+ preloss-1]) == nothing
        #     continue
        # end
        wage_path[ind,1:last_index-first_index+1]    = income[set[ind][1], first_index:last_index]
    end
    return wage_path
end

function consumption_loss(sim::Simulations,preloss::Int64=5,postloss::Int64=8)
    @unpack T = Primitives()
    @unpack emp_streak, c = sim
    matrix = emp_streak[:,preloss+1:end]  #omitting the first period where everyone is unemployed
    set = findall(iszero,matrix)
    c_path = zeros(size(set,1),preloss+postloss+1)
    for ind in 1:size(set,1) 
        first_index         = set[ind][2]
        last_index          = min(set[ind][2] + preloss + postloss,T)
        # if income[set[ind][1], set[ind][2]+ preloss-1] == 0
        #     continue
        # end
        # if findall(isone,wage[set[ind][1], 1:set[ind][2]+ preloss-1]) == nothing
        #     continue
        # end
        c_path[ind,1:last_index-first_index+1]    = c[set[ind][1], first_index:last_index]
    end
    return c_path
end