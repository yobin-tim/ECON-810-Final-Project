# Load packages
using Parameters, Statistics, Distributions, ProgressBars, SharedArrays, Distributed

# Include strutures with primitives and results
include("structures.jl")
addprocs(4) 
# Initialize the model
@everywhere function Init()
    
    # Initialize the primitives
    #! Everithing that should be calibrated I will pass as a parameter to the primitive struct
    # A                   = 1.2                       #productivity scale #todo: if possible, calibrate                 
    # μ                   = 1 - 0.2                   #todo: calibrate if possible
    # prim                = Primitives(A,μ)
    prim                = Primitives()
    # Initialize the value functions
    U                   = SharedArray{Float64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    W_L                 = SharedArray{Float64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    W_H                 = SharedArray{Float64}(prim.n_kPoints, prim.n_hPoints, prim.n_SPoints, prim.T)
    #Pre-compute last period value function
    for khS in 1:(prim.n_kPoints * prim.n_hPoints * prim.n_SPoints) 
        k,h,S           = Tuple(CartesianIndices((prim.n_kPoints, prim.n_hPoints, prim.n_SPoints))[khS])
        U[k,h,S,prim.T]      = prim.u.(prim.k_grid[k] .* (1 + prim.r) .+ prim.b)
        W_L[k,h,S,prim.T]    = prim.u.(prim.k_grid[k] .* (1 + prim.r) .+ prim.R_L(prim.T) .* prim.h_grid[h])
        W_H[k,h,S,prim.T]    = prim.u.(prim.k_grid[k] .* (1 + prim.r) .+ prim.R_H(prim.T) .* prim.h_grid[h])
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
    K_net = (1 + prim.r) * repeat(prim.k_grid', prim.n_kPoints,1) .- repeat(prim.k_grid, 1,prim.n_kPoints) # K_net[i,j]=(1+r)*k[i]-k'[j]
    ℓ_effective = repeat(prim.h_grid', prim.n_sPoints,1) .* (1 .- repeat(prim.s_grid, 1,prim.n_hPoints))# ℓ_effective[i,j]=h[i]*(1-s[j])
    hs = repeat(prim.h_grid, 1, prim.n_sPoints) + ( repeat(prim.h_grid, 1 , prim.n_sPoints) .* repeat(prim.s_grid', prim.n_hPoints, 1)).^prim.α # hs[i,j]=h[i]+(h[i]*s[j])^α
    h_next = zeros(prim.n_hPoints, prim.n_sPoints, prim.n_zPoints) # h_next[i,j,k]=exp(z[k])*h[i]+(h[i]*s[j])^α
    # Multiply random shock to hs
    for z_i in 1:prim.n_zPoints
        h_next[:,:,z_i] = exp(prim.z_grid[z_i]) .* hs
    end
    h_next[ h_next .<  prim.h_min ] .= prim.h_min # Replace everything below h_min with h_min
    h_next[ h_next .>  prim.h_max ] .= prim.h_max # Replace everything above h_max with h_max
    h_next = round.(h_next, digits = 1) # Round to nearest integer
    h_next_indexes = zeros(size(h_next))
    for i in 1:prim.n_hPoints
        for j in 1:prim.n_sPoints
            for k in 1:prim.n_zPoints
                # h_next_indexes[i,j,k] = find(prim.h_grid == h_next[i,j,k])
                h_next_indexes[i,j,k]     = findfirst(x -> x == h_next[i,j,k], prim.h_grid)[1]
            end
        end
    end

    pre_comp            = Pre_Computed(K_net', ℓ_effective', h_next_indexes)
    # Return primitives and results
    return prim, res, pre_comp
end

# Compute the value function
function vfn(prim::Primitives, res::Results, pre_comp::Pre_Computed)
    
    # Unpack primitives, results and pre-computed variables
    @unpack T, n_kPoints, n_hPoints, n_SPoints, R_L, R_H, h_grid, k_grid, r, u, β, s_grid, S_grid, n_sPoints, z_trProb, z_grid, H, h_min, h_max, α, b, Π, S_bar, δ, μ  = prim
    @unpack val_fun, k_pol, s_pol = res
    @unpack K_net, ℓ_effective, h_next_indexes = pre_comp
    # Compute the value function

    for t in ProgressBar(T-1:-1:1) 
        @sync @distributed for khS in 1:(n_kPoints * n_hPoints * n_SPoints)
            k,h,S       = Tuple(CartesianIndices((n_kPoints,n_hPoints,n_SPoints))[khS])
            # Skip unnatainable states
            if S > (t-1)*(n_sPoints-1)+1
                continue
            end
            Sprime_ind = S:S+n_sPoints-1 # index of the next period's S for any given s
            h_next_ind = h_next_indexes[h,:,:] # (fixed h) index of the next period's h for any given s all shocks z
            income_WL   = R_L(t) .* ℓ_effective[h, :] # Possible income for each selection of s (Firm L)
            income_WH   = R_H(t) .* ℓ_effective[h, :] # Possible income for each selection of s (Firm H)
            k_net       = K_net[k,:] # Net capital for each possible choice of k_next
            
            C_WL   = repeat(k_net, 1 ,n_sPoints) .+ repeat(income_WL', n_kPoints, 1) # Consumption for each selection of (k', s) (Firm L)
            C_WH   = repeat(k_net, 1 ,n_sPoints) .+ repeat(income_WH', n_kPoints, 1) # Consumption for each selection of (k', s) (Firm H)
            C_U = b .+ repeat(k_net, 1 ,n_sPoints)
            U_WL = zeros(size(C_WL)) # Utility for each selection of (k', s) (Firm L)
            U_WH = zeros(size(C_WH)) # Utility for each selection of (k', s) (Firm H)
            U_U = zeros(size(C_U))
            U_WL[C_WL .< 0] .= -Inf # Set impossible consumption to -Inf
            U_WL[C_WL .>= 0] .= u.(C_WL[C_WL .>= 0]) # Set possible consumption to utility
            U_WH[C_WH .< 0] .= -Inf # Set impossible consumption to -Inf
            U_WH[C_WH .>= 0] .= u.(C_WH[C_WH .>= 0]) # Set possible consumption to utility
            U_U[C_U .< 0] .= -Inf # Set impossible consumption to -Inf
            U_U[C_U .>= 0] .= u.(C_U[C_U .>= 0]) # Set possible consumption to utility
            
            # cand_val_U  = -Inf
            # cand_val_WL = -Inf
            # cand_val_WH = -Inf

            temp_exp_WL = zeros(n_kPoints, n_sPoints)
            temp_exp_WH = zeros(n_kPoints, n_sPoints)
            temp_exp_U  = zeros(n_kPoints, n_sPoints)
            for kprime in 1:n_kPoints
                for s in 1:n_sPoints
                    temp_exp_WL[kprime,s] = (val_fun.W_L[kprime, h_next_ind, Sprime_ind[s], t+1] * z_trProb)[1]
                    temp_exp_WH[kprime,s] = (val_fun.W_H[kprime, h_next_ind, Sprime_ind[s], t+1] * z_trProb)[1]
                    temp_exp_U[kprime,s]  = (val_fun.U[kprime, h_next_ind, Sprime_ind[s], t+1] * z_trProb  )[1]
                end
            end
            # println(size(temp_exp_WL))
            # println(size(U_WL))
            temp_WL = U_WL + β*( (1-δ) * temp_exp_WL + δ * temp_exp_U )
            temp_WH = U_WH + β*( (1-δ) * temp_exp_WH + δ * temp_exp_U ) 
            Π_sh = Π.(1 .- s_grid, S_grid[S], t)'
            if S_grid[S] < S_bar
                temp_U = U_U .+ β*( μ * Π_sh .* temp_exp_WL .+ (1 .- μ * Π_sh) .* temp_exp_U)
            else
                temp_U = U_U .+ β*( Π_sh .* ( μ .* temp_exp_WL + (1 - μ) .* temp_exp_WH ) .+ (1 .- Π_sh) .* temp_exp_U)
            end

            # Update value functions 
            val_fun.W_L[k,h,S,t] = maximum(temp_WL)
            val_fun.W_H[k,h,S,t] = maximum(temp_WH)
            val_fun.U[k,h,S,t]   = maximum(temp_U)

            # Update policy functions
            k_pol.W_L[k,h,S,t], s_pol.W_L[k,h,S,t] = Tuple(argmax(temp_WL))
            k_pol.W_H[k,h,S,t], s_pol.W_H[k,h,S,t] = Tuple(argmax(temp_WH))
            k_pol.U[k,h,S,t],   s_pol.U[k,h,S,t]   = Tuple(argmax(temp_U))

        end # end of khS loop
    end # End of t loop 
end # End of function


@time prim, res, pre_comp = Init();
@time vfn(prim, res, pre_comp)

#             for s in 1:(n_sPoints)
#                 hprime = h_next[h,s,:] 
#                 indexes         = zeros(Int64,size(hprime,1),1)
#                 for i in 1:size(hprime,1)
#                     indexes[i]     = findfirst(x -> x == hprime[i], h_grid)[1]
#                 end
#                 Sprime          = round.(min(S_grid[S] + s_grid[s], S_grid[end]), digits = 1)
#                 Sprime_ind      =  findfirst(x -> x == Sprime, S_grid)[1]
#                 #Loop for U
#                 for kprime in 1:n_kPoints
#                     c_U = budget_U .- k_grid[kprime]
#                     if c_U < 0
#                         break
#                     end
#                     if S_grid[S] < S_bar
#                         val_func_U = u(c_U) .+ β * z_trProb' *  ( Π(1 .-s_grid[s], S_grid[Sprime_ind],t) .* W_L[kprime,indexes,Sprime_ind,t+1] .+ (1 .- Π(1 .-s_grid[s], S_grid[Sprime_ind],t) .* U[kprime, indexes, Sprime_ind, t+1]))
#                     else
#                         val_func_U = u(c_U) .+ β * z_trProb' *  ( Π(1 .-s_grid[s], S_grid[Sprime_ind],t) .* ( μ .* W_L[kprime,indexes,Sprime_ind,t+1] .+ (1 - μ) .* W_H[kprime,indexes,Sprime_ind,t+1] ) .+ (1 .- Π(1 .-s_grid[s], S_grid[Sprime_ind],t) .* U[kprime, indexes, Sprime_ind, t+1]))
#                     end
#                     if val_func_U[1] > cand_val_U
#                         U[k,h,S,t]              = val_func_U[1]
#                         k_pol_U[k,h,S,t]        = k_grid[kprime]
#                         s_pol_U[k,h,S,t]        = s_grid[s]
#                         k_pol_ind_U[k,h,S,t]    = kprime
#                         s_pol_ind_U[k,h,S,t]    = s
#                         cand_val_U              = val_func_U[1]
#                     end
#                 end
#                 #Loop for W_L and W_H
#                 for kprime in 1:n_kPoints
#                     c_WL = budget_WL[s] .- k_grid[kprime]
#                     c_WH = budget_WH[s] .- k_grid[kprime]
#                     if c_WH < 0 # since c_WL is always lower than c_WH
#                         break
#                     end
#                     val_func_WL = u(c_WL) .+ β * z_trProb' * ( (1 .- δ) .* W_L[kprime, indexes, Sprime_ind, t+1] .+ δ .* U[kprime, indexes, Sprime_ind, t+1] )
#                     val_func_WH = u(c_WL) .+ β * z_trProb' * ( (1 .- δ) .* W_H[kprime, indexes, Sprime_ind, t+1] .+ δ .* U[kprime, indexes, Sprime_ind, t+1] )
#                     if val_func_WL[1] > cand_val_WL
#                         W_L[k,h,S,t]                = val_func_WL[1]
#                         k_pol_W_L[k,h,S,t]          = k_grid[kprime]
#                         s_pol_W_L[k,h,S,t]          = s_grid[s]
#                         k_pol_ind_W_L[k,h,S,t]      = kprime
#                         s_pol_ind_W_L[k,h,S,t]      = s
#                         cand_val_WL                 = val_func_WL[1]
#                     end
#                     if val_func_WH[1] > cand_val_WH
#                         W_H[k,h,S,t]                = val_func_WH[1]
#                         k_pol_W_H[k,h,S,t]          = k_grid[kprime]
#                         s_pol_W_H[k,h,S,t]          = s_grid[s]
#                         k_pol_ind_W_H[k,h,S,t]      = kprime
#                         s_pol_ind_W_H[k,h,S,t]      = s
#                         cand_val_WH                 = val_func_WH[1]
#                     end
#                 end
#             end
#         end
#     end
# end