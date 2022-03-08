# Load packages
@everywhere using Parameters, Statistics, Distributions, ProgressBars, SharedArrays

# Include strutures with primitives and results
include("structures.jl")

# Initialize the model
@everywhere function Init()
    
    # Initialize the primitives
    #! Everithing that should be calibrated I will pass as a parameter to the primitive struct
    A                   = 1.2                       #productivity scale #todo: if possible, calibrate                 
    μ                   = 1 - 0.2                   #todo: calibrate if possible
    prim                = Primitives(A,μ)
    
    # Initilize the value functions
    U                   = SharedArray{Float64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    W_L                 = SharedArray{Float64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    W_H                 = SharedArray{Float64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    #Pre-compute last period value function
    for khS in 1:(prim.n_kPoints * prim.n_hPoints * prim.n_SPoints) 
        k,h,S           = Tuple(CartesianIndices((n_kPoints,n_hPoints,n_SPoints))[khS])
        U[k,h,S,T]      = prim.u.(prim.k_grid[k] .* (1 + prim.r) .+ prim.b)
        W_L[k,h,S,T]    = prim.u.(prim.k_grid[k] .* (1 + prim.r) .+ prim.R_L(prim.T) .* prim.h_grid[h])
        W_H[k,h,S,T]    = prim.u.(prim.k_grid[k] .* (1 + prim.r) .+ prim.R_H(prim.T) .* prim.h_grid[h])
    end
    val_fun             = ValueFunction(U, W_L, W_H)
    
    # Initialize the asset holdings policy functions
    k_pol_U             = SharedArray{Float64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    k_pol_ind_U         = SharedArray{Int64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    k_pol_W_L           = SharedArray{Float64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    k_pol_ind_W_L       = SharedArray{Int64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    k_pol_W_H           = SharedArray{Float64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    k_pol_ind_W_H       = SharedArray{Int64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    k_pol               = PolicyFunctionAssets(k_pol_U, k_pol_ind_U, k_pol_W_L, k_pol_ind_W_L, k_pol_W_H, k_pol_ind_W_H)
    
    # Initialize the schooling choice policy functions
    s_pol_U             = SharedArray{Float64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    s_pol_ind_U         = SharedArray{Int64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    s_pol_W_L           = SharedArray{Float64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    s_pol_ind_W_L       = SharedArray{Int64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    s_pol_W_H           = SharedArray{Float64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    s_pol_ind_W_H       = SharedArray{Int64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    s_pol               = PolicyFunctionSchooling(s_pol_U, s_pol_ind_U, s_pol_W_L, s_pol_ind_W_L, s_pol_W_H, s_pol_ind_W_H)
    
    # Initialize the structure to hold the results
    res                 = Results(val_fun, k_pol, s_pol)
    
    # Pre-compute auxiliar variables
    K_net = (1 + prim.r) * repeat(prim.k_grid, prim.n_kPoints,1) .- repeat(prim.k_grid', 1,prim.n_kPoints) # K_net[i,j]=(1+r)*k[i]-k'[j]
    ℓ_effective = repeat(prim.h_grid, prim.n_hPoints,1) .* (1 .- repeat(prim.s_grid', 1,prim.n_sPoints))# ℓ_effective[i,j]=h[i]*(1-s[j])
    hs = repeat(prim.h_grid, prim.n_hPoints,1) + ( repeat(prim.h_grid, prim.n_hPoints,1) .* repeat(prim.s_grid', 1,prim.n_sPoints)).^prim.α # hs[i,j]=h[i]+(h[i]*s[j])^α
    h_next = zeros(prim.n_kPoints, prim.n_hPoints, prim.n_zPoints) # h_next[i,j,k]=exp(z[k])*h[i]+(h[i]*s[j])^α
    # Multiply random shock to hs
    for z_i in 1:prim.n_zPoints
        h_next[:,:,z_i] = exp(prim.z_grid[z_i]) .* hs
    end
    h_next[ h_next .<  prim.h_min ] .= prim.h_min # Replace everything below h_min with h_min
    h_next[ h_next .>  prim.h_max ] .= prim.h_max # Replace everything above h_max with h_max
    h_next = round.(h_next, digits = 1) # Round to nearest integer
    pre_comp            = Pre_Computed(K_net, ℓ_effective, h_next)
    # Return primitives and results
    return prim, res, pre_comp
end

# Compute the value function
function vfn(prim::Primitives, res::Results, pre_comp::Pre_Computed)
    
    # Unpack primitives, results and pre-computed variables
    @unpack T, n_kPoints, n_hPoints, n_SPoints, R_L, R_H, h_grid, k_grid, r, u, β, s_grid, S_grid, n_sPoints, z_trProb, z_grid, H, h_min, h_max, α, b, Π, S_bar, δ, μ  = prim
    @unpack val_fun, k_pol, s_pol = res
    @unpack K_net, ℓ_effective, h_next = pre_comp

    for t in ProgressBar(T-1:-1:1) 
        @sync @distributed for khS in 1:(n_kPoints * n_hPoints * n_SPoints)
            k,h,S       = Tuple(CartesianIndices((n_kPoints,n_hPoints,n_SPoints))[khS])
            if S > t*n_sPoints
                continue
            end
            income_WL   = R_L(t) .* ℓ_effective[h, :]
            income_WH   = R_H(t) .* ℓ_effective[h, :]
            # budget_U    = b .+ k_grid[k] * (1 + r)
            cand_val_U  = -Inf
            cand_val_WL = -Inf
            cand_val_WH = -Inf
            for s in 1:(n_sPoints)
                hprime = h_next[h,s,:]
                # hprime          = round.(exp.(z_grid) .* H(h_grid[h],s_grid[s]), digits = 1)  #use the law of motion
                # hprime          = replace(x -> x > h_max ? h_max : x, hprime)   #replace values exceeding upper bound by h_max
                # hprime          = replace(x -> x < h_min ? h_min : x, hprime)   #replace values below lower bound by h_min 
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