function tauchenMethod(μ::Float64, σ::Float64, ρ::Float64, z_num::Int64; q::Int64=3, tauchenoptions=0, dshift=0) 
    
    # Tauchen's method for generating Markov chains. 
    # This is based on TauchenMethod.m from the VFIToolkit-matlab toolbox by Robert Kirkby
    # https://github.com/vfitoolkit/VFIToolkit-matlab

    # Create z_vals vector and transition matrix for the discrete markov process approximation of AR(1) process z'=μ+ρ*z+e, e~N(0,σ²), by Tauchens method
    # Inputs
    #   μ              - AR(1) process z'=μ+ρ*z+ε, ε~N(0,σ²)
    #   ρ              - AR(1) process z'=μ+ρ*z+ε, ε~N(0,σ²)
    #   σ²              - AR(1) process z'=μ+ρ*z+ε, ε~N(0,σ²)
    #   z_num           - number of z_vals in discretization of z (must be an odd number)
    #   q              - max number of std devs from mean (default=3)
    # Optional Inputs
    # Outputs
    #   z_vals         - column vector containing the z_num states of the discrete approximation of z
    #   Π    - transition matrix of the discrete approximation of z;
    #                    Π(i,j) is the probability of transitioning from state i to state j
    #
    # Helpful info: Π
    #   Var(z)=σ²/(1-ρ^2); note that if μ=0, then σ²z=σ²/(1-ρ^2).
    ###############
    
    
    if z_num==1
        z_vals=[μ/(1-ρ)]; #expected value of z
        Π=[1];

        return z_vals, Π
    end

        # σ = sqrt(σ²); #stddev of ε
        z_star = μ/(1-ρ) #expected value of z
        σ_z = σ/sqrt(1-ρ^2) #stddev of z
        z = z_star*ones(z_num, 1) + range(-q*σ_z, stop = q*σ_z, length = z_num)  
        ω = z[2] - z[1] #Note that all the points are equidistant by construction.
        
        zi=z*ones(1,z_num);

        zj=dshift*ones(z_num,z_num)+ones(z_num,1)*z'
        
        dist = Normal(μ, σ)

        P_part1 = cdf( dist, zj .+ ω/2-ρ * zi)
        P_part2 = cdf( dist, zj .- ω/2-ρ * zi)
        
        P = P_part1 - P_part2
        P[:, 1] = P_part1[:, 1]
        P[:, z_num] = 1 .- P_part2[:,z_num]
        

        z_vals=z;
        Π=P; #(z,zprime)
        return z_vals, Π
end