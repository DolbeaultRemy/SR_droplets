include("Functions_time_evol_boxes.jl")

# Nbr of Threads
println("Nbr of Threads = $(Threads.nthreads())")

### Define the system ###
# Nbr of particles
N_list = [2:2:10;]
r = 100 # Nbr of repetitions

# Normalisation parameters
λ = 421e-9
γ = 32.7e6 # In Hz

# Physical values
ω0 = (2π*ustrip(c_0)/λ)
ωl = ω0
kl = [ustrip(c_0)/ωl, 0, 0] # Laser along x
Ω_RO = 1e7 # Taken from Barbut arXiv:2412.02541v1

# Fixed density
n0 = 1000 #2e3 # atoms per unit of volume (already normalized)
d0_lb = 5e-10 # Minimum distance between the atoms (lower boundary) in m

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
window_t, window_var = 2, 4
threshold_box = 1e-5; # Time over which the boxes are computed

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
popup_t_N, list_t_N, nbr_error_t_N, sol_t_N = [], [], [], []

for (i, N) in ProgressBar(enumerate(N_list))
    @load "op_list/op_list_$N.jdl2" op_list
    tasks = map(chunks) do chunk # Split the different distributions into chuncks solved on each core
        Threads.@spawn solve_random_distrib(chunk, functs, op_list, N, n0, d0_lb, window_t, window_var, threshold_box, i, false)
    end
    sol_tasks = fetch.(tasks)
    push!(list_t_N, vcat([s[1] for s in sol_tasks]...))
    push!(popup_t_N, vcat([s[2] for s in sol_tasks]...))
    push!(nbr_error_t_N, vcat([s[3] for s in sol_tasks]...))
    push!(sol_t_N, vcat([s[4] for s in sol_tasks]...))
end

@save "solutions/sol_N_$(N_list)_r_$(r).jld2" list_t_N popup_t_N nbr_error_t_N sol_t_N