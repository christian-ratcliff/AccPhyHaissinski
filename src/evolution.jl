"""
evolution.jl - Particle evolution and beam dynamics

This file implements the core longitudinal evolution algorithm and particle generation.
It handles the multi-turn tracking of particles through the accelerator,
including RF cavity effects, synchrotron radiation, and collective effects.
"""

using StructArrays
using StaticArrays
using Statistics
using Random
using Distributions
using ProgressMeter
using LoopVectorization

"""
    generate_particles(
        μ_z::T, μ_E::T, σ_z::T, σ_E::T, num_particles::Int,
        energy::T, mass::T, ϕs::T, freq_rf::T
    ) where T<:Float64 -> Tuple{StructArray{Particle{T}}, T, T, T}

Generate initial particle distribution.
"""



function generate_particles(
    μ_z::T, μ_E::T, σ_z::T, σ_E::T, num_particles::Int,
    energy::T, mass::T, ϕs::T, freq_rf::T
    ) where {T<:Union{Real, StochasticAD.StochasticTriple}}

    # Initial sampling
    initial_sample_size = min(10000, num_particles)
    z_samples = rand(Normal(μ_z, σ_z), initial_sample_size)
    E_samples = rand(Normal(μ_E, σ_E), initial_sample_size)

    # Covariance matrix
    Σ = Symmetric([cov(z_samples, z_samples) cov(z_samples, E_samples);
                   cov(z_samples, E_samples) cov(E_samples, E_samples)])

    μ_vec = SVector{2,eltype(μ_z)}(μ_z, μ_E)
    dist_total = MvNormal(μ_vec, Σ)

    # Relativistic factors
    γ = energy / mass
    β = sqrt(1 - 1/γ^2)
    rf_factor = freq_rf * 2π / (β * SPEED_LIGHT)

    # Sample
    samples = rand(dist_total, num_particles)
    z_vals = samples[1, :]
    ΔE_vals = samples[2, :]

    # Promote to StochasticTriple if needed
    if T <: StochasticAD.StochasticTriple
        z_vals = StochasticAD.stochastic_triple.(z_vals)
        ΔE_vals = StochasticAD.stochastic_triple.(ΔE_vals)
    end

    # Create particles
    coords = Coordinate.(z_vals, ΔE_vals)
    particles = StructArray{Particle{typeof(coords[1].z)}}((StructArray(coords),))

    return particles, σ_E, σ_z, energy
end

function longitudinal_evolve!(
    particles_float64::StructArray{Particle{Float64}},
    E0_param::TE, 
    mass_param::TM, 
    voltage_param::TV, 
    harmonic_param::TH,
    radius_param::TR, 
    α_c_param::TA, 
    ϕs_param::TP, 
    freq_rf_param::TF,
    n_turns::Int
) where {TE, TM, TV, TH, TR, TA, TP, TF}
    
    # Determine the working type from physics parameters
    sample_computation = E0_param + voltage_param + radius_param + α_c_param + ϕs_param + freq_rf_param
    T = typeof(sample_computation)
    
    # Create zero of the working type T through arithmetic
    zero_T = sample_computation * 0
    
    # Copy Float64 particles to working type T using arithmetic operations
    n_particles = length(particles_float64)
    coords_T = Vector{Coordinate{T}}(undef, n_particles)
    
    for i in 1:n_particles
        # Use arithmetic to create T-typed values, avoiding type conversion
        z_T = zero_T + particles_float64.coordinates.z[i]
        ΔE_T = zero_T + particles_float64.coordinates.ΔE[i]
        coords_T[i] = Coordinate{T}(z_T, ΔE_T)
    end
    
    # Create StructArray with working type T
    particles = StructArray{Particle{T}}((StructArray(coords_T),))
    
    # Pre-compute physical constants in working type T
    γ0 = E0_param / mass_param
    β0 = sqrt(1 - 1/γ0^2)
    η0 = α_c_param - 1/γ0^2
    sin_ϕs = sin(ϕs_param)
    rf_factor = freq_rf_param * 2π / (β0 * SPEED_LIGHT)
    
    # Initial spreads - convert to working type through arithmetic
    σ_E0_initial = zero_T + std(particles_float64.coordinates.ΔE)
    
    # Main evolution loop - now all arrays can hold type T
    for turn in 1:n_turns
        # RF voltage kick
        for i in 1:n_particles
            ϕ_val = -particles.coordinates.z[i] * rf_factor + ϕs_param
            particles.coordinates.ΔE[i] += voltage_param * (sin(ϕ_val) - sin_ϕs)
        end
        
        # Quantum excitation (stochastic effect)
        ∂U_∂E = 4 * 8.85e-5 * (E0_param/1e9)^3 / radius_param
        excitation = sqrt(1-(1-∂U_∂E)^2) * σ_E0_initial
        for i in 1:n_particles
            # randn() gives Float64, arithmetic with excitation promotes to T
            particles.coordinates.ΔE[i] += excitation * randn()
        end
        
        # Synchrotron radiation damping
        damping_factor = 1 - ∂U_∂E
        for i in 1:n_particles
            particles.coordinates.ΔE[i] *= damping_factor
        end
        
        # Update reference energy
        E0_param = E0_param - 4 * 8.85e-5 * (E0_param/1e9)^3 / radius_param * E0_param / 4
        γ0 = E0_param / mass_param
        β0 = sqrt(1 - 1/γ0^2)
        η0 = α_c_param - 1/γ0^2
        
        # Update phase advance
        coeff = 2π * harmonic_param * η0 / (β0^2 * E0_param)
        for i in 1:n_particles
            ϕ_i = -(particles.coordinates.z[i] * rf_factor - ϕs_param)
            ϕ_i += coeff * particles.coordinates.ΔE[i]
            particles.coordinates.z[i] = (-ϕ_i + ϕs_param) / rf_factor
        end
        
        # Update RF factor
        rf_factor = freq_rf_param * 2π / (β0 * SPEED_LIGHT)
    end
    
    # Return final energy spread in MeV
    σ_E_final = std(particles.coordinates.ΔE)
    return σ_E_final / 1e6
end


function longitudinal_evolve_enzyme_wrapper(
    particles_float64::StructArray{Particle{Float64}},
    param_value::Float64,
    param_type::Symbol,
    fixed_params::NamedTuple,
    n_turns::Int,
    rng_seed::Int
)
    # Set the random seed for reproducible stochastic effects
    Random.seed!(rng_seed)
    
    # Unpack fixed parameters and replace the one we're differentiating
    if param_type == :voltage
        return longitudinal_evolve_parametric(
            particles_float64, fixed_params.E0_ini, fixed_params.mass, 
            param_value, fixed_params.harmonic, fixed_params.radius, 
            fixed_params.α_c, fixed_params.ϕs, fixed_params.freq_rf, n_turns
        )
    elseif param_type == :energy
        return longitudinal_evolve_parametric(
            particles_float64, param_value, fixed_params.mass, 
            fixed_params.voltage, fixed_params.harmonic, fixed_params.radius, 
            fixed_params.α_c, fixed_params.ϕs, fixed_params.freq_rf, n_turns
        )
    else
        error("Unknown parameter type: $param_type")
    end
end


