# Packages for plotting and saving/loading data 
using NPZ, Plots

include("functions.jl")

# Plot style
begin
    theme(:vibrant); #:vibrant
    default(fontfamily="Computer Modern", framestyle=:box); # LaTex-style
    gr(size = (400, 400)); # default plot size
end

# This function to save data to files
function save_data(res, path, suffix)
    path = path * suffix * "/";
    run(`mkdir -p $path`)
    npzwrite(path*"U"*suffix*".npz", res.val_fun.U)
    npzwrite(path*"W_L"*suffix*".npz", res.val_fun.W_L)
    npzwrite(path*"W_H"*suffix*".npz", res.val_fun.W_H)
    npzwrite(path*"k_pol_U"*suffix*".npz", res.k_pol.U)
    npzwrite(path*"s_pol_U"*suffix*".npz", res.s_pol.U)
    npzwrite(path*"k_pol_WL"*suffix*".npz", res.k_pol.W_L)
    npzwrite(path*"s_pol_WL"*suffix*".npz", res.s_pol.W_L)
    npzwrite(path*"k_pol_WH"*suffix*".npz", res.k_pol.W_H)
    npzwrite(path*"s_pol_WH"*suffix*".npz", res.s_pol.W_H)
end

# This function read data from files and returns all the structures needed

function read_data(path, suffix)

    # First initialize the structures
    prim, res, pre_comp = Init()
    path = path * suffix * "/";
    res.val_fun.U = npzread(path*"U"*suffix*".npz")
    res.val_fun.W_L = npzread(path*"W_L"*suffix*".npz")
    res.val_fun.W_H = npzread(path*"W_H"*suffix*".npz")
    res.k_pol.U = npzread(path*"k_pol_U"*suffix*".npz"  )
    res.s_pol.U = npzread(path*"s_pol_U"*suffix*".npz"  )
    res.k_pol.W_L = npzread(path*"k_pol_WL"*suffix*".npz")
    res.s_pol.W_L = npzread(path*"s_pol_WL"*suffix*".npz")
    res.k_pol.W_H = npzread(path*"k_pol_WH"*suffix*".npz")
    res.s_pol.W_H = npzread(path*"s_pol_WH"*suffix*".npz")

    return prim, res, pre_comp
end