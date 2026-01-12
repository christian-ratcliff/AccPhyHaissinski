function setup_physics_parameters()
    """Set up base physical parameters for the simulation"""
    E0_ini = 4e9
    mass = MASS_ELECTRON
    voltage = 5e6
    harmonic = 360
    radius = 250.0
    pipe_radius = 0.00025
    α_c = 3.68e-4
    ϕs = 5π/6
    freq_rf = let
        γ = E0_ini/mass
        β = sqrt(1 - 1/γ^2)
        (ϕs + 10*π/180) * β * SPEED_LIGHT / (2π)
    end

    # Distribution parameters
    μ_z = 0.0
    μ_E = 0.0
    σ_E0 = 1e6
    σ_z0 = let
        γ = E0_ini/mass
        β = sqrt(1 - 1/γ^2)
        ω_rev = 2 * π / ((2*π*radius) / (β*SPEED_LIGHT))
        sqrt(2 * π) * SPEED_LIGHT / ω_rev * sqrt(α_c*E0_ini/harmonic/voltage/abs(cos(ϕs))) * σ_E0 / E0_ini
    end

    return (E0_ini=E0_ini, mass=mass, voltage=voltage, harmonic=harmonic, 
            radius=radius, pipe_radius=pipe_radius, α_c=α_c, ϕs=ϕs, freq_rf=freq_rf,
            μ_z=μ_z, μ_E=μ_E, σ_E0=σ_E0, σ_z0=σ_z0)
end