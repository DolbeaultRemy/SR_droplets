include("Functions_time_evol_boxes_density_stripe.jl")

# Nbr of Threads
println("Nbr of Threads = $(Threads.nthreads())")

### Define the system ###
# Nbr of particles
N = 10
sat = [40.1;]
# vcat([0.1:0.1:0.9;], [1:35;]) # [1:35;] #[0.1:5:35.1;]
pathname_density_distrib = "Groundstate_Bext_90_deg_as_98.h5"
n_list = [10:50:1000;]  # List of densities to simulate.
save_sol = false # Save the full time evolution or just the SS
r = 100 # Nbr of repetitions


# Normalisation parameters
λ = 421e-9
γ = 32.7e6 # In Hz

# Physical values
ω0 = (2π*ustrip(c_0)/λ)
ωl = ω0
kl = [2π/λ, 0, 0] # Laser along x
Ω_RO = γ*sqrt.(sat/2)

# Normalization
ω0 = ω0 / γ
ωl = ωl / γ
kl = kl * λ 
Ω_RO = Ω_RO / γ



# Quantization axis along z
e = [0, 0, 1.]

# Integration parameter
t_end = 1e3
window_t, window_var = 2, 5 # Times over which the boxes are computed. The first one is to smooth popup, the second one is the length of the boxes used to find the SS
threshold_box = 5e-4; # Threshold for the boxes

# Load the set of differential equations
fsolve_diff_eq(du, u, params) = ccall(("diffeqf", "libs/liballfuncs_$N.dll"), Cvoid, (Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}), du, u, params)

# Load operators list
@load "op_list/op_list_$N.jdl2" op_list

Create_dir(sat, n_list)

list_r = 1:r
chunks = Iterators.partition(list_r, cld(length(list_r), Threads.nthreads()))

for (i, n) in ProgressBar(enumerate(n_list))
    tasks = map(chunks) do chunk # Split the different repetitions into chuncks solved on each core
        Threads.@spawn solve_random_distrib_vary_Ω(chunk, fsolve_diff_eq, op_list, N, n_list, window_t, window_var, threshold_box, t_end, save_sol, Ω_RO[1])
    end
    sol_tasks = fetch.(tasks)
    @save "solutions/Sat_$(minimum(sat))to$(maximum(sat))_n0_$(round(mean(n_list), digits=2))_N_$N/sol_N_$(N)_Sat_$(sat[1])_stripe_distribution_$(pathname_density_distrib[1:end-3])_n0_$(n).jld2" N n_list sat sol_tasks
    sol_tasks = nothing
end