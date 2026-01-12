function run_enzyme_ad_benchmark(params, param_type::Symbol, x, enzyme_samples, n_particles, n_turns)
    """Run Enzyme AD benchmark with stochastic sampling"""
    
    derivatives = Float64[]
    
    # Create fixed parameters for enzyme wrapper
    fixed_params = if param_type == :voltage
        (E0_ini=params.E0_ini, mass=params.mass, harmonic=params.harmonic,
         radius=params.radius, α_c=params.α_c, ϕs=params.ϕs, freq_rf=params.freq_rf,
         voltage=params.voltage)
    elseif param_type == :energy
        (mass=params.mass, voltage=params.voltage, harmonic=params.harmonic,
         radius=params.radius, α_c=params.α_c, ϕs=params.ϕs, freq_rf=params.freq_rf,
         E0_ini=params.E0_ini)
    else
        error("Unknown parameter type: $param_type")
    end
    
    # Get the base parameter value to scale
    base_value = param_type == :voltage ? params.voltage : params.E0_ini
    param_value = x * base_value
    
    # Time the Enzyme AD calculation
    enzyme_time = @elapsed begin
        for i in 1:enzyme_samples
            # Generate fresh particles for each sample to capture stochastic variation
            particles, _, _, _ = generate_particles(
                params.μ_z, params.μ_E, params.σ_z0, params.σ_E0, 
                n_particles, params.E0_ini, params.mass, params.ϕs, params.freq_rf
            )
            
            # Use unique seed for each sample to ensure different stochastic behavior
            rng_seed = 1000 + i
            
            # Create wrapper function for this specific run
            wrapper_func(p) = longitudinal_evolve_enzyme_wrapper(
                particles, p, param_type, fixed_params, n_turns, rng_seed
            )
            
            # Use Enzyme forward mode autodiff
            deriv = Enzyme.autodiff(Enzyme.Forward, Const(wrapper_func), Enzyme.Duplicated(param_value, 1.0))[1]
            
            # Convert to derivative w.r.t. scale factor (chain rule)
            # d/d(scale) = d/d(param) * d(param)/d(scale) = deriv * base_value
            scale_deriv = deriv * base_value
            
            push!(derivatives, scale_deriv)
        end
    end
    
    enzyme_mean = mean(derivatives)
    enzyme_err = std(derivatives) / sqrt(enzyme_samples)
    
    return enzyme_mean, enzyme_err, enzyme_time, enzyme_samples  # 1 function eval per sample
end