include("Functions_time_evol_boxes_density_stripe.jl")

# Nbr of Threads
println("Nbr of Threads = $(Threads.nthreads())")

### Define the system ###
# Nbr of particles
N = 10
sat = [0.1:5:35.1;]
pathname_density_distrib = "Groundstate_Bext_90_deg_as_98.h5"
@load "Density_distributions\\$(pathname_density_distrib[1:end-3]).jld2" CS_densities_cuts_interpolated_norm
n_list = CS_densities_cuts_interpolated_norm#[1:1:length(CS_densities_cuts_interpolated_norm)] # List of density to simulate
save_sol = false # Save the full time evolution or just the SS

# Normalisation parameters
λ = 421e-9
γ = 32.7e6 # In Hz

# Physical values
ω0 = (2π*ustrip(c_0)/λ)
ωl = ω0
kl = [2π/λ, 0, 0] # Laser along x
Ω_RO = γ*sqrt.(sat/2)

d0_lb = 1e-10 # Minimum distance between the atoms (lower boundary) in m

# Normalization
ω0 = ω0 / γ
ωl = ωl / γ
kl = kl * λ 
Ω_RO = Ω_RO / γ
d0_lb = d0_lb / λ

# Quantization axis along z
e = [0, 0, 1.]

# Integration parameter
t_end = 1e3
window_t, window_var = 2, 5 # Times over which the boxes are computed. The first one is to smooth popup, the second one is the length of the boxes used to find the SS
threshold_box = 5e-4; # Threshold for the boxes

# Load the set of differential equations
fsolve_diff_eq(du, u, params) = ccall(("diffeqf", "libs/liballfuncs_$N.dll"), Cvoid, (Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}), du, u, params)

### Computation ###
list_n = 1:length(n_list)
chunks_n = Iterators.partition(list_n, cld(length(list_n), Threads.nthreads()))

# Load operators list
@load "op_list/op_list_$N.jdl2" op_list

for (i, Ω) in ProgressBar(enumerate(Ω_RO))
    tasks = map(chunks_n) do chunk # Split the different density distributions into chuncks solved on each core
        Threads.@spawn solve_random_distrib_vary_Ω(chunk, fsolve_diff_eq, op_list, N, n_list, d0_lb, window_t, window_var, threshold_box, t_end, save_sol, Ω)
    end
    sol_tasks = fetch.(tasks)
    @save "solutions/sol_N_$(N)_sat_$(sat[i])_stripe_distribution_$(pathname_density_distrib).jld2" sol_tasks
end