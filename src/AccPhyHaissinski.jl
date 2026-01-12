module AccPhyHaissinski


include("constants.jl")
include("data_structures.jl")
include("utils.jl")
include("rf_kick.jl")
include("synchrotron_radiation.jl")
include("evolution.jl")
include("quantum_excitation.jl")
include("params.jl")
include("wakefield.jl")

include("benchmarks/evolution_enzyme.jl")

# Export constants

export SPEED_LIGHT, ELECTRON_CHARGE, MASS_ELECTRON

# Export data structures
export Coordinate, Particle, BeamTurn, SimulationParameters, SimulationBuffers

# Export core functions
export generate_particles, longitudinal_evolve!
export quantum_excitation!, synchrotron_radiation!, apply_wakefield_inplace!, rf_kick!, synchrotron_radiation!

export run_enzyme_ad_benchmark

# Export utilities
export z_to_ϕ, ϕ_to_z, calc_rf_factor, create_simulation_buffers,copyto_particles!, next_power_of_two

end
