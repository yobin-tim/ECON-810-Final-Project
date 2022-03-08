@everywhere using Parameters, Statistics, Distributions, SharedArrays
include("tauchen.jl")

@everywhere @with_kw struct Primitives
    β::Float64                  = 0.99
    r::Float64                  = 0.04
    T::Int64                    = 30                        #in years
    A::Float64                  = 1.2                    #productivity scale                  
    R_L::Function               = (t) -> (1.0019)^(t-1)     #rental rate on labor
    R_H::Function               = (t) -> A * (1.0019)^(t-1)
    μ_z::Float64                = -0.029
    σ_z::Float64                = sqrt(0.11)
    σ::Float64                  = 2
    z                           = Normal(μ_z, σ_z)
    n_zPoints::Int64            = 8
    tauchen                     = tauchenMethod(μ_z,σ_z,0.0,n_zPoints)
    z_grid                      = tauchen[1]
    z_trProb                    = tauchen[2][1,:]
    b::Float64                  = 0.5
    δ::Float64                  = 0.1
    μ::Float64                  = 1 - 0.2                     #todo: calibrate if possible
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
    threshold::Float64          = 0.4
    S_bar::Float64              = S_grid[round(Int64,n_SPoints * threshold)]                           #todo: calibrate, if possible
    #Simulations
    nSim::Int64                 = 1000
end

# Structure to hold the value functions
@everywhere mutable struct ValueFunction
    U                   ::SharedArray{Float64,4}        # ValueFunction for an unemployed individual
    W_L                 ::SharedArray{Float64,4}        # Value function for an employed individual in a Low-type firm
    W_H                 ::SharedArray{Float64,4}        # Value function for an employed individual in a High-type firm
end

# Structures to hold the policy functions
# Asset holdings policy functions
@everywhere mutable struct PolicyFunctionAssets
    U             ::SharedArray{Float64,4}        # Asset holdings for individuals that are uemployed
    ind_U         ::SharedArray{Int64,  4}        # (index) Asset holdings for individuals that are uemployed
    W_L           ::SharedArray{Float64,4}        # Asset holdings for individuals that are employed in the Low-type firm
    ind_W_L       ::SharedArray{Int64,  4}        # (index) Asset holdings for individuals that are employed in the Low-type firm
    W_H           ::SharedArray{Float64,4}        # Asset holdings for individuals that are employed in the High-type firm
    ind_W_H       ::SharedArray{Int64,  4}        # (index) Asset holdings for individuals that are employed in the High-type firm
end 



# Schooling choice policy functions
@everywhere mutable struct PolicyFunctionSchooling
    U             ::SharedArray{Float64,4}        # Schooling choice for individuals that are uemployed
    ind_U         ::SharedArray{Int64,  4}        # (index) Schooling choice for individuals that are uemployed
    W_L           ::SharedArray{Float64,4}        # Schooling choice for individuals that are employed in the Low-type firm
    ind_W_L       ::SharedArray{Int64,  4}        # (index)Schooling choice for individuals that are employed in the Low-type firm
    W_H           ::SharedArray{Float64,4}        # Schooling choice for individuals that are employed in the High-type firm
    ind_W_H       ::SharedArray{Int64,  4}        # (index) Schooling choice for individuals that are employed in the High-type firm
end

# Structures to hold the resutls
@everywhere mutable struct Results
    val_fun             ::ValueFunction                 # Value functions
    k_pol               ::PolicyFunctionAssets          # Asset holdings policy functions
    s_pol               ::PolicyFunctionSchooling       # Schooling choice policy functions
end

# Auxiliar structures to hold pre-computed values
# @everywhere mutable struct Consumption
#     C_U                 ::SharedArray{Float64,2}        # Consumption unemployed individuals
#     C_W_L               ::SharedArray{Float64,5}        # Consumption employed in the Low-type firm
#     C_W_H               ::SharedArray{Float64,2}        # Consumption employed in the High-type firm
# end

@everywhere mutable struct Pre_Computed
    K_net                 ::Array{Float64,2}        # Net of assets: (1+r)k - k'
    ℓ_effective           ::Array{Float64,2}        # Effective labor: h(1-s) 
    h_next_indexes        ::Array{Int64,3}        # Next period human capital
end