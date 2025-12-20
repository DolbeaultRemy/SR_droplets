include("./Functions_time_evol_decay.jl")

# Nbr of Threads
println("Nbr of Threads = $(Threads.nthreads())")

### Define the system ###
# Nbr of particles
N_list = [4]#[2:2:10;]
r = 5#10 # Nbr of repetitions

# Normalisation parameters
λ = 421e-9
γ = 32.7e6 # In Hz

# Physical values
ω0 = (2π*ustrip(c_0)/λ)
ωl = ω0

# Fixed density
n0 = 1000 #2e3 # atoms per unit of volume (already normalized)
d0_lb = 5e-10 # Minimum distance between the atoms (lower boundary) in m

# Normalization
ω0 = ω0 / γ
ωl = ωl / γ

# Quantization axis along z
e = [0, 0, 1.]

# Integration parameter
t_end = 1e1

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
    @load "../SS/op_list/op_list_$N.jdl2" op_list
    tasks = map(chunks) do chunk # Split the different distributions into chuncks solved on each core
        Threads.@spawn solve_decay_SS(chunk, functs, N, r, i, true)
    end
    sol_tasks = fetch.(tasks)
    @save "solutions/sol_N_$(N)_r_$(r).jld2" sol_tasks
end