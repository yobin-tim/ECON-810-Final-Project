@everywhere using Parameters, Statistics, Distributions, ProgressBars, SharedArrays

include("../Problem Set 2/tauchen.jl")

@everywhere @with_kw struct Primitives
    β::Float64                  = 0.99
    r::Float64                  = 0.04
    T::Int64                    = 30                        #in years
    R::Function                 = (t) -> (1.0019)^(t-1)     #rental rate on labor
    μ_z::Float64                = -0.029
    σ_z::Float64                = sqrt(0.11)
    σ::Float64                  = 2
    z                           = Normal(μ_z, σ_z)
    n_zPoints::Int64            = 10
    tauchen                     = tauchenMethod(μ_z,σ_z,0.0,n_zPoints)
    z_grid                      = tauchen[1]
    z_trProb                    = tauchen[2][1,:]
    α::Float64                  = 0.70
    H::Function                 = (h,s) -> h + (h*s)^α
    u::Function                 = (c) -> (c^(1 - σ)-1)/(1 - σ)
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
    k_max::Float64              = 60.0
    n_kPoints::Int64            = 121
    k_grid::Array{Float64,1}    = range(k_min, k_max, length=n_kPoints)
    #Share of time to invest in human capital
    n_sPoints::Int64            = 21
    s_grid::Array{Float64,1}    = range(0.0,1.0,length = n_sPoints)
    #Simulations
    nSim::Int64                 = 1000
end

@everywhere mutable struct Results
    V::SharedArray{Float64,3}
    k_pol::SharedArray{Float64,3}
    k_pol_ind::SharedArray{Int64,3}
    s_pol::SharedArray{Float64,3}
    s_pol_ind::SharedArray{Int64,3}
end


@everywhere function Init()
    prim        = Primitives()
    V           = SharedArray{Float64}(prim.n_kPoints, prim.n_hPoints, prim.T)
    k_pol       = SharedArray{Float64}(prim.n_kPoints, prim.n_hPoints, prim.T)
    k_pol_ind   = SharedArray{Int64}(prim.n_kPoints, prim.n_hPoints, prim.T)
    s_pol       = SharedArray{Float64}(prim.n_kPoints, prim.n_hPoints, prim.T)
    s_pol_ind   = SharedArray{Int64}(prim.n_kPoints, prim.n_hPoints, prim.T)
    res         = Results(V,k_pol, k_pol_ind, s_pol, s_pol_ind)
    return prim, res
end

@everywhere function tr_mat(h::Float64, s::Float64)
    @unpack α, H, z_grid, z_trProb, h_grid, h_max = Primitives()
    hprime          = round.(exp.(z_grid) .* H(h,s), digits = 1)  #use the law of motion
    hprime          = replace(x -> x > h_max ? h_max : x, hprime)   #replace values exceeding upper bound by  
    h_trProb        = zeros(size(hprime,1),1)
    indexes         = zeros(Int64,size(hprime,1),1)
    for i in 1:size(hprime,1)
        indexes[i]     = findfirst(x -> x == hprime[i], h_grid)[1]
        h_trProb[i]    = z_trProb[i]
    end
    return hprime, h_trProb, indexes #h_prime grid may not be required
end

function vfn(prim::Primitives, res::Results)
    @unpack T, n_kPoints, n_hPoints, R, h_grid, k_grid, r, u, β, s_grid, n_sPoints  = prim
    @unpack V, k_pol, s_pol, k_pol_ind, s_pol_ind = res

    V[:,:,T] = u.(k_grid .* (1 + r) .+ R(T) .* h_grid')
    for t in ProgressBar(T-1:-1:1) 
        @sync @distributed for kh in 1:(n_kPoints * n_hPoints)
            k,h = Tuple(CartesianIndices((n_kPoints,n_hPoints))[kh])
            # for h in 1:n_hPoints
                budget = R(t) .* h_grid[h] .*  (1 .- s_grid) .+ k_grid[k]  * (1 + r)
                cand_val = -Inf
                @sync @distributed for skprime in 1:(n_kPoints *n_sPoints)
                    kprime, s = Tuple(CartesianIndices((n_kPoints, n_sPoints))[skprime]) 
                    c = budget[s] - k_grid[kprime]
                    if c < 0
                        break
                    end
                    hprime_grid, h_trProb, hprime_ind = tr_mat(h_grid[h],s_grid[s])
                    val_func = u(c) .+ β * h_trProb' *  V[kprime,hprime_ind,t+1]
                    if val_func[1] > cand_val
                        V[k,h,t]            = val_func[1]
                        k_pol[k,h,t]        = k_grid[kprime]
                        k_pol_ind[k,h,t]    = kprime
                        s_pol[k,h,t]        = s_grid[s]
                        s_pol_ind[k,h,t]    = s
                    end
                end
            # end
        end
    end
end

#vfn2 not using tr_mat and reducing computation of h' given (h,s).
function vfn2(prim::Primitives, res::Results)
    @unpack T, n_kPoints, n_hPoints, R, h_grid, k_grid, r, u, β, s_grid, n_sPoints, z_trProb, z_grid, H, h_min, h_max, α  = prim
    @unpack V, k_pol, s_pol, k_pol_ind, s_pol_ind = res

    V[:,:,T] = u.(k_grid .* (1 + r) .+ R(T) .* h_grid')
    for t in ProgressBar(T-1:-1:1) 
        @sync @distributed for kh in 1:(n_kPoints * n_hPoints)
            k,h = Tuple(CartesianIndices((n_kPoints,n_hPoints))[kh])
            # for h in 1:n_hPoints
                budget = R(t) .* h_grid[h] .*  (1 .- s_grid) .+ k_grid[k]  * (1 + r)
                cand_val = -Inf
                for s in 1:(n_sPoints)
                    hprime          = round.(exp.(z_grid) .* H(h_grid[h],s_grid[s]), digits = 1)  #use the law of motion
                    hprime          = replace(x -> x > h_max ? h_max : x, hprime)   #replace values exceeding upper bound by h_max
                    hprime          = replace(x -> x < h_min ? h_min : x, hprime)   #replace values below lower bound by h_min 
                    h_trProb        = zeros(size(hprime,1),1)
                    indexes         = zeros(Int64,size(hprime,1),1)
                    for i in 1:size(hprime,1)
                        indexes[i]     = findfirst(x -> x == hprime[i], h_grid)[1]
                        h_trProb[i]    = z_trProb[i]
                    end
                    for kprime in 1:n_kPoints
                        c = budget[s] - k_grid[kprime]
                        if c < 0
                            break
                        end
                        
                        val_func = u(c) .+ β * h_trProb' *  V[kprime,indexes,t+1]
                        if val_func[1] > cand_val
                            V[k,h,t]            = val_func[1]
                            k_pol[k,h,t]        = k_grid[kprime]
                            s_pol[k,h,t]        = s_grid[s]
                            k_pol_ind[k,h,t]    = kprime
                            s_pol_ind[k,h,t]    = s
                        end
                    end
                    
                end
            # end
        end
    end
end

# for kprime in 1:n_kPoints, s in 1:n_sPoints 
#     c = budget[s] - k_grid[k]
#     if c < 0
#         break
#     end
#     hprime_grid, h_trProb, hprime_ind = tr_mat(h_grid[h],s_grid[s])
#     val_func = u(c) .+ β * h_trProb' *  V[k,hprime_ind,t+1]
#     if val_func[1] > cand_val
#         V[k,h,t]        = val_func[1]
#         k_pol[k,h,t]    = k_grid[kprime]
#         s_pol[k,h,t]    = s_grid[s]
#     end
# end

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

t1 = plot(1:4, [2 3 4 5.25]',
        xlims  = (1,4),
        xticks = 1:1:4,
        ylims  = (0,6),
        legend = false,
        xlabel = "Time to Maturity (years)",
        ylabel = "Yield to Maturity (%)"
        )