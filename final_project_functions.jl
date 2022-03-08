@everywhere using Parameters, Statistics, Distributions, ProgressBars, SharedArrays

include("tauchen.jl")

@everywhere @with_kw struct Primitives
    β::Float64                  = 0.99
    r::Float64                  = 0.04
    T::Int64                    = 30                        #in years
    A::Float64                  = 1.2                       #productivity scale #todo: if possible, calibrate                 
    R_L::Function               = (t) -> (1.0019)^(t-1)     #rental rate on labor
    R_H::Function               = (t) -> (A * 1.0019)^(t-1)
    μ_z::Float64                = -0.029
    σ_z::Float64                = sqrt(0.11)
    σ::Float64                  = 2
    z                           = Normal(μ_z, σ_z)
    n_zPoints::Int64            = 8
    tauchen                     = tauchenMethod(μ_z,σ_z,0.0,n_zPoints)
    z_grid                      = tauchen[1]
    z_trProb                    = tauchen[2][1,:]
    b::Float64                  = 0.2
    δ::Float64                  = 0.1
    μ::Float64                  = 0.2                   #todo: calibrate if possible
    α::Float64                  = 0.70
    H::Function                 = (h,s) -> h + (h*s)^α
    u::Function                 = (c) -> (c^(1 - σ)-1)/(1 - σ)
    Π::Function                 = (γ, S, t) -> γ * (S/t)
    #human capital grid
    h_min::Float64              = 0.1
    h_max::Float64              = 10.0
    n_hPoints::Int64            = 100
    h_grid::Array{Float64,1}    = range(h_min, h_max, length=n_hPoints)
    #Initial human capital
    μ_hc::Float64               = 2.0
    σ_hc::Float64               = sqrt(0.5)
    hc_i                        = Normal(μ_hc, σ_hc)
    #Asset grid 
    k_min::Float64              = 0.0
    k_max::Float64              = 40.0
    n_kPoints::Int64            = 41
    k_grid::Array{Float64,1}    = range(k_min, k_max, length=n_kPoints)
    #Share of time to invest in human capital
    n_sPoints::Int64            = 6
    s_grid::Array{Float64,1}    = range(0.0,1.0,length = n_sPoints)
    #Cumulative years of schooling
    n_SPoints::Int64            = (n_sPoints-1) * T + 1
    S_grid::Array{Float64,1}    = range(0.0, 1.0 * T, length = n_SPoints)
    S_bar::Float64              = S_grid[121]                           #todo: calibrate, if possible
    #Simulations
    nSim::Int64                 = 1000
end

@everywhere mutable struct Results
    U::SharedArray{Float64,4}
    W_L::SharedArray{Float64,4}
    W_H::SharedArray{Float64,4}
    k_pol_U::SharedArray{Float64,4}
    k_pol_ind_U::SharedArray{Int64,4}
    s_pol_U::SharedArray{Float64,4}
    s_pol_ind_U::SharedArray{Int64,4}
    k_pol_W_L::SharedArray{Float64,4}
    k_pol_ind_W_L::SharedArray{Int64,4}
    s_pol_W_L::SharedArray{Float64,4}
    s_pol_ind_W_L::SharedArray{Int64,4}
    k_pol_W_H::SharedArray{Float64,4}
    k_pol_ind_W_H::SharedArray{Int64,4}
    s_pol_W_H::SharedArray{Float64,4}
    s_pol_ind_W_H::SharedArray{Int64,4}
end


@everywhere function Init()
    prim        = Primitives()
    U           = SharedArray{Float64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    W_L         = SharedArray{Float64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    W_H         = SharedArray{Float64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    k_pol_U       = SharedArray{Float64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    k_pol_ind_U   = SharedArray{Int64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    s_pol_U       = SharedArray{Float64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    s_pol_ind_U   = SharedArray{Int64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    k_pol_W_L       = SharedArray{Float64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    k_pol_ind_W_L   = SharedArray{Int64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    s_pol_W_L       = SharedArray{Float64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    s_pol_ind_W_L   = SharedArray{Int64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    k_pol_W_H       = SharedArray{Float64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    k_pol_ind_W_H   = SharedArray{Int64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    s_pol_W_H       = SharedArray{Float64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    s_pol_ind_W_H   = SharedArray{Int64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    res         = Results(U, W_L, W_H ,k_pol_U, k_pol_ind_U, s_pol_U, s_pol_ind_U,k_pol_W_L, k_pol_ind_W_L, s_pol_W_L, s_pol_ind_W_L,k_pol_W_H, k_pol_ind_W_H, s_pol_W_H, s_pol_ind_W_H)
    return prim, res
end


function vfn2(prim::Primitives, res::Results)
    @unpack T, n_kPoints, n_hPoints, n_SPoints, R_L, R_H, h_grid, k_grid, r, u, β, s_grid, S_grid, n_sPoints, z_trProb, z_grid, H, h_min, h_max, α, b, Π, S_bar, δ, μ  = prim
    @unpack U, W_L, W_H ,k_pol_U, k_pol_ind_U, s_pol_U, s_pol_ind_U,k_pol_W_L, k_pol_ind_W_L, s_pol_W_L, s_pol_ind_W_L,k_pol_W_H, k_pol_ind_W_H, s_pol_W_H, s_pol_ind_W_H = res

    #Update last period value function
    for khS in 1:(n_kPoints * n_hPoints * n_SPoints) 
        k,h,S = Tuple(CartesianIndices((n_kPoints,n_hPoints,n_SPoints))[khS])
        U[k,h,S,T]      = u.(k_grid[k] .* (1 + r) .+ b)
        W_L[k,h,S,T]    = u.(k_grid[k] .* (1 + r) .+ R_L(T) .* h_grid[h])
        W_H[k,h,S,T]    = u.(k_grid[k] .* (1 + r) .+ R_H(T) .* h_grid[h])
    end
    

    for t in ProgressBar(T-1:-1:1) 
        @sync @distributed for khS in 1:(n_kPoints * n_hPoints * n_SPoints)
            k,h,S       = Tuple(CartesianIndices((n_kPoints,n_hPoints,n_SPoints))[khS])
            budget_WL   = R_L(t) .* h_grid[h] .*  (1 .- s_grid) .+ k_grid[k]  * (1 + r)
            budget_WH   = R_H(t) .* h_grid[h] .*  (1 .- s_grid) .+ k_grid[k]  * (1 + r)
            budget_U    = b .+ k_grid[k] * (1 + r)
            cand_val_U  = -Inf
            cand_val_WL = -Inf
            cand_val_WH = -Inf
            if S_grid[S] < S_bar
                for s in 1:(n_sPoints)
                    hprime          = round.(exp.(z_grid) .* H(h_grid[h],s_grid[s]), digits = 1)  #use the law of motion
                    hprime          = replace(x -> x > h_max ? h_max : x, hprime)   #replace values exceeding upper bound by h_max
                    hprime          = replace(x -> x < h_min ? h_min : x, hprime)   #replace values below lower bound by h_min 
                    indexes         = zeros(Int64,size(hprime,1),1)
                    for i in 1:size(hprime,1)
                        indexes[i]     = findfirst(x -> x == hprime[i], h_grid)[1]
                    end
                    Sprime          = round.(min(S_grid[S] + s_grid[s], S_grid[end]), digits = 1)
                    Sprime_ind      =  findfirst(x -> x == Sprime, S_grid)[1]
                    #Loop for U
                    for kprime in 1:n_kPoints
                        c_U = budget_U .- k_grid[kprime]
                        if c_U < 0
                            break
                        end
                        val_func_U = u(c_U) .+ β * z_trProb' *  ( Π(1 .-s_grid[s], S_grid[Sprime_ind],t) .* W_L[kprime,indexes,Sprime_ind,t+1] .+ (1 .- Π(1 .-s_grid[s], S_grid[Sprime_ind],t) .* U[kprime, indexes, Sprime_ind, t+1]))
                        if val_func_U[1] > cand_val_U
                            U[k,h,S,t]              = val_func_U[1]
                            k_pol_U[k,h,S,t]        = k_grid[kprime]
                            s_pol_U[k,h,S,t]        = s_grid[s]
                            k_pol_ind_U[k,h,S,t]    = kprime
                            s_pol_ind_U[k,h,S,t]    = s
                            cand_val_U              = val_func_U[1]
                        end
                    end
                    #Loop for W_L and W_H
                    for kprime in 1:n_kPoints
                        c_WL = budget_WL[s] .- k_grid[kprime]
                        c_WH = budget_WH[s] .- k_grid[kprime]
                        if c_WH < 0 # since c_WL is always lower than c_WH
                            break
                        end
                        val_func_WL = u(c_WL) .+ β * z_trProb' * ( (1 .- δ) .* W_L[kprime, indexes, Sprime_ind, t+1] .+ δ .* U[kprime, indexes, Sprime_ind, t+1] )
                        val_func_WH = u(c_WL) .+ β * z_trProb' * ( (1 .- δ) .* W_H[kprime, indexes, Sprime_ind, t+1] .+ δ .* U[kprime, indexes, Sprime_ind, t+1] )
                        if val_func_WL[1] > cand_val_WL
                            W_L[k,h,S,t]                = val_func_WL[1]
                            k_pol_W_L[k,h,S,t]          = k_grid[kprime]
                            s_pol_W_L[k,h,S,t]          = s_grid[s]
                            k_pol_ind_W_L[k,h,S,t]      = kprime
                            s_pol_ind_W_L[k,h,S,t]      = s
                            cand_val_WL                 = val_func_WL[1]
                        end
                        if val_func_WH[1] > cand_val_WH
                            W_H[k,h,S,t]                = val_func_WH[1]
                            k_pol_W_H[k,h,S,t]          = k_grid[kprime]
                            s_pol_W_H[k,h,S,t]          = s_grid[s]
                            k_pol_ind_W_H[k,h,S,t]      = kprime
                            s_pol_ind_W_H[k,h,S,t]      = s
                            cand_val_WH                 = val_func_WH[1]
                        end
                    end
                end
            else
                for s in 1:(n_sPoints)
                    hprime          = round.(exp.(z_grid) .* H(h_grid[h],s_grid[s]), digits = 1)  #use the law of motion
                    hprime          = replace(x -> x > h_max ? h_max : x, hprime)   #replace values exceeding upper bound by h_max
                    hprime          = replace(x -> x < h_min ? h_min : x, hprime)   #replace values below lower bound by h_min 
                    indexes         = zeros(Int64,size(hprime,1),1)
                    for i in 1:size(hprime,1)
                        indexes[i]     = findfirst(x -> x == hprime[i], h_grid)[1]
                    end
                    Sprime          = round.(min(S_grid[S] + s_grid[s], S_grid[end]), digits = 1)
                    Sprime_ind      =  findfirst(x -> x == Sprime, S_grid)[1]
                    #Loop for U
                    for kprime in 1:n_kPoints
                        c_U = budget_U[s] .- k_grid[kprime]
                        if c_U < 0
                            break
                        end
                        val_func_U = u(c_U) .+ β * z_trProb' *  ( Π(1 .-s_grid[s], S_grid[Sprime_ind],t) .* ( μ .* W_L[kprime,indexes,Sprime_ind,t+1] .+ (1 - μ) .* W_H[kprime,indexes,Sprime_ind,t+1] ) .+ (1 .- Π(1 .-s_grid[s], S_grid[Sprime_ind],t) .* U[kprime, indexes, Sprime_ind, t+1]))
                        if val_func_U[1] > cand_val_U
                            U[k,h,S,t]              = val_func_U[1]
                            k_pol_U[k,h,S,t]        = k_grid[kprime]
                            s_pol_U[k,h,S,t]        = s_grid[s]
                            k_pol_ind_U[k,h,S,t]    = kprime
                            s_pol_ind_U[k,h,S,t]    = s
                            cand_val_U              = val_func_U[1]
                        end
                    end
                    #Loop for W_L and W_H
                    for kprime in 1:n_kPoints
                        c_WL = budget_WL[s] - k_grid[kprime]
                        c_WH = budget_WH[s] - k_grid[kprime]
                        if c_WH < 0 # since c_WL is always lower than c_WH
                            break
                        end
                        val_func_WL = u(c_WL) .+ β * z_trProb' * ( (1 .- δ) .* W_L[kprime, indexes, Sprime_ind, t+1] .+ δ .* U[kprime, indexes, Sprime_ind, t+1] )
                        val_func_WH = u(c_WL) .+ β * z_trProb' * ( (1 .- δ) .* W_H[kprime, indexes, Sprime_ind, t+1] .+ δ .* U[kprime, indexes, Sprime_ind, t+1] )
                        if val_func_WL[1] > cand_val_WL
                            W_L[k,h,S,t]                = val_func_WL[1]
                            k_pol_W_L[k,h,S,t]          = k_grid[kprime]
                            s_pol_W_L[k,h,S,t]          = s_grid[s]
                            k_pol_ind_W_L[k,h,S,t]      = kprime
                            s_pol_ind_W_L[k,h,S,t]      = s
                            cand_val_WL                 = val_func_WL[1]
                        end
                        if val_func_WH[1] > cand_val_WH
                            W_H[k,h,S,t]                = val_func_WH[1]
                            k_pol_W_H[k,h,S,t]          = k_grid[kprime]
                            s_pol_W_H[k,h,S,t]          = s_grid[s]
                            k_pol_ind_W_H[k,h,S,t]      = kprime
                            s_pol_ind_W_H[k,h,S,t]      = s
                            cand_val_WH                 = val_func_WH[1]
                        end
                    end
                end
            end
        end
    end
end

function vfn(prim::Primitives, res::Results)
    @unpack T, n_kPoints, n_hPoints, n_SPoints, R_L, R_H, h_grid, k_grid, r, u, β, s_grid, S_grid, n_sPoints, z_trProb, z_grid, H, h_min, h_max, α, b, Π, S_bar, δ, μ  = prim
    @unpack U, W_L, W_H ,k_pol_U, k_pol_ind_U, s_pol_U, s_pol_ind_U,k_pol_W_L, k_pol_ind_W_L, s_pol_W_L, s_pol_ind_W_L,k_pol_W_H, k_pol_ind_W_H, s_pol_W_H, s_pol_ind_W_H = res

    #Update last period value function
    for khS in 1:(n_kPoints * n_hPoints * n_SPoints) 
        k,h,S = Tuple(CartesianIndices((n_kPoints,n_hPoints,n_SPoints))[khS])
        U[k,h,S,T]      = u.(k_grid[k] .* (1 + r) .+ b)
        W_L[k,h,S,T]    = u.(k_grid[k] .* (1 + r) .+ R_L(T) .* h_grid[h])
        W_H[k,h,S,T]    = u.(k_grid[k] .* (1 + r) .+ R_H(T) .* h_grid[h])
    end
    

    for t in ProgressBar(T-1:-1:1) 
        @sync @distributed for khS in 1:(n_kPoints * n_hPoints * n_SPoints)
            k,h,S       = Tuple(CartesianIndices((n_kPoints,n_hPoints,n_SPoints))[khS])
            budget_WL   = R_L(t) .* h_grid[h] .*  (1 .- s_grid) .+ k_grid[k]  * (1 + r)
            budget_WH   = R_H(t) .* h_grid[h] .*  (1 .- s_grid) .+ k_grid[k]  * (1 + r)
            budget_U    = b .+ k_grid[k] * (1 + r)
            cand_val_U  = -Inf
            cand_val_WL = -Inf
            cand_val_WH = -Inf
            for s in 1:(n_sPoints)
                hprime          = round.(exp.(z_grid) .* H(h_grid[h],s_grid[s]), digits = 1)  #use the law of motion
                hprime          = replace(x -> x > h_max ? h_max : x, hprime)   #replace values exceeding upper bound by h_max
                hprime          = replace(x -> x < h_min ? h_min : x, hprime)   #replace values below lower bound by h_min 
                indexes         = zeros(Int64,size(hprime,1),1)
                for i in 1:size(hprime,1)
                    indexes[i]     = findfirst(x -> x == hprime[i], h_grid)[1]
                end
                Sprime          = round.(min(S_grid[S] + s_grid[s], S_grid[end]), digits = 1)
                Sprime_ind      =  findfirst(x -> x == Sprime, S_grid)[1]
                #Loop for U
                for kprime in 1:n_kPoints
                    c_U = budget_U .- k_grid[kprime]
                    if c_U < 0
                        break
                    end
                    if S_grid[S] < S_bar
                        val_func_U = u(c_U) .+ β * z_trProb' *  ( Π(1 .-s_grid[s], S_grid[Sprime_ind],t) .* W_L[kprime,indexes,Sprime_ind,t+1] .+ (1 .- Π(1 .-s_grid[s], S_grid[Sprime_ind],t) .* U[kprime, indexes, Sprime_ind, t+1]))
                    else
                        val_func_U = u(c_U) .+ β * z_trProb' *  ( Π(1 .-s_grid[s], S_grid[Sprime_ind],t) .* ( μ .* W_L[kprime,indexes,Sprime_ind,t+1] .+ (1 - μ) .* W_H[kprime,indexes,Sprime_ind,t+1] ) .+ (1 .- Π(1 .-s_grid[s], S_grid[Sprime_ind],t) .* U[kprime, indexes, Sprime_ind, t+1]))
                    end
                    if val_func_U[1] > cand_val_U
                        U[k,h,S,t]              = val_func_U[1]
                        k_pol_U[k,h,S,t]        = k_grid[kprime]
                        s_pol_U[k,h,S,t]        = s_grid[s]
                        k_pol_ind_U[k,h,S,t]    = kprime
                        s_pol_ind_U[k,h,S,t]    = s
                        cand_val_U              = val_func_U[1]
                    end
                end
                #Loop for W_L and W_H
                for kprime in 1:n_kPoints
                    c_WL = budget_WL[s] .- k_grid[kprime]
                    c_WH = budget_WH[s] .- k_grid[kprime]
                    if c_WH < 0 # since c_WL is always lower than c_WH
                        break
                    end
                    val_func_WL = u(c_WL) .+ β * z_trProb' * ( (1 .- δ) .* W_L[kprime, indexes, Sprime_ind, t+1] .+ δ .* U[kprime, indexes, Sprime_ind, t+1] )
                    val_func_WH = u(c_WL) .+ β * z_trProb' * ( (1 .- δ) .* W_H[kprime, indexes, Sprime_ind, t+1] .+ δ .* U[kprime, indexes, Sprime_ind, t+1] )
                    if val_func_WL[1] > cand_val_WL
                        W_L[k,h,S,t]                = val_func_WL[1]
                        k_pol_W_L[k,h,S,t]          = k_grid[kprime]
                        s_pol_W_L[k,h,S,t]          = s_grid[s]
                        k_pol_ind_W_L[k,h,S,t]      = kprime
                        s_pol_ind_W_L[k,h,S,t]      = s
                        cand_val_WL                 = val_func_WL[1]
                    end
                    if val_func_WH[1] > cand_val_WH
                        W_H[k,h,S,t]                = val_func_WH[1]
                        k_pol_W_H[k,h,S,t]          = k_grid[kprime]
                        s_pol_W_H[k,h,S,t]          = s_grid[s]
                        k_pol_ind_W_H[k,h,S,t]      = kprime
                        s_pol_ind_W_H[k,h,S,t]      = s
                        cand_val_WH                 = val_func_WH[1]
                    end
                end
            end
        end
    end
end

@everywhere mutable struct Simulations
    hc::Array{Float64,2}
    hc_ind::Array{Int64,2}
    s_inv::Array{Float64,2}
    s_ind::Array{Int64,2}
    k::Array{Float64,2}
    k_ind::Array{Int64,2}
    c::Array{Float64,2}
    wage_income::Array{Float64,2}
end


function Init2(prim::Primitives)
    hc          = zeros(prim.nSim, prim.T)
    hc_ind      = ones(prim.nSim, prim.T)
    s_inv       = zeros(prim.nSim, prim.T)
    s_ind       = ones(prim.nSim, prim.T)
    k           = zeros(prim.nSim, prim.T)
    k_ind       = ones(prim.nSim, prim.T)
    c           = zeros(prim.nSim, prim.T)
    wage_income      = zeros(prim.nSim, prim.T)
    sim = Simulations(hc, hc_ind, s_inv, s_ind, k, k_ind, c, wage_income)
    return sim
end


function runsim(prim::Primitives, res::Results, sim::Simulations)
    @unpack hc, hc_ind, s_inv, s_ind, k, k_ind, c, wage_income = sim
    @unpack V, k_pol, s_pol, k_pol_ind, s_pol_ind = res
    @unpack T, nSim, hc_i, n_kPoints, n_hPoints, R, h_grid, k_grid, r, u, β, s_grid, n_sPoints, z_trProb, z_grid, H, h_min, h_max, α  = prim
    hc_ind[:,1] = rand(1:n_hPoints, nSim)                   #initial human capital index, round off to match grid
    hc[:,1]     = h_grid[hc_ind[:,1]]               
    for t in 1:T 
        for i in 1:nSim
            s_inv[i,t]          = s_pol[k_ind[i,t],hc_ind[i,t],t]
            s_ind[i,t]          = s_pol_ind[k_ind[i,t],hc_ind[i,t],t]
            wage_income[i,t]    = R(t) * hc[i,t] * (1 - s_inv[i])
            budget              = wage_income[i,t] + k[i,t] * (1 + r)
            if t < T
                k[i,t+1]            = k_pol[k_ind[i,t],hc_ind[i,t],t]
                k_ind[i,t+1]        = k_pol_ind[k_ind[i,t],hc_ind[i,t],t]
                c[i,t]              = budget - k[i,t+1]
                hc[i,t+1]           = round.(rand(z_grid) * H(hc[i,t],s_inv[i,t]), digits = 1)
                if hc[i,t+1] > h_max
                    hc[i,t+1] = h_max
                elseif hc[i,t+1] < h_min
                    hc[i,t+1] = h_min
                end 
                hc_ind[i,t+1]       = findfirst(x -> x == hc[i,t+1], h_grid)[1]
            else
                c[i,t]              = budget
            end
            
            
        end
    end
end