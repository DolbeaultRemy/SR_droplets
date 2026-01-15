include("Functions_time_evol_boxes.jl")

# Nbr of Threads
println("Nbr of Threads = $(Threads.nthreads())")

### Define the system ###
# Nbr of particles
N_list = [10;]
r = 100 # Nbr of repetitions

# Normalisation parameters
λ = 421e-9
γ = 32.7e6 # In Hz

# Physical values
ω0 = (2π*ustrip(c_0)/λ)
ωl = ω0
kl = [2π/λ, 0, 0] # Laser along x
# sat = [0:10:100;]
sat = [0:1e-2:1;]
# sat = [5;]
Ω_RO = γ*sqrt.(sat/2)

# Fixed density
n0 = 1e2 #2e3 # atoms per unit of volume (already normalized)
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
t_end = 1e2
window_t, window_var = 2, 4 # Times over which the boxes are computed. The first one is to smooth popup, the second one is the length of the boxes used to find the SS
threshold_box = 1e-3; # Threshold for the boxes

# Prepare the wrapper
const N_FUNCS = length(N_list)  # Total function nbr
const functs = Vector{Function}(undef, N_FUNCS)

for (i, N) in enumerate(N_list)
    libpath ="libs/liballfuncs_$N.dll"
    functs[i] = load_f("diffeqf", libpath)
end



### Computation ###
list_r = 1:r
chunks = Iterators.partition(list_r, cld(length(list_r), Threads.nthreads()))

for (i, N) in ProgressBar(enumerate(N_list))
    @load "op_list/op_list_$N.jdl2" op_list
    for (j, Ω) in ProgressBar(enumerate(Ω_RO))
        tasks = map(chunks) do chunk # Split the different distributions into chuncks solved on each core
            Threads.@spawn solve_random_distrib_vary_Ω(chunk, functs, op_list, N, n0, d0_lb, window_t, window_var, threshold_box, i, true, Ω)
        end
        sol_tasks = fetch.(tasks)
        @save "solutions/sol_N_$(N)_r_$(r)_sat_$(sat[j])_n0_$n0.jld2" sol_tasks
    end
end