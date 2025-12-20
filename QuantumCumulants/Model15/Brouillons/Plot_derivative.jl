include("MiniCollectiveSpins.jl")
include("Functions_time_evol_decay_rate.jl")
using PyPlot
using Statistics
using JLD2

# Nbr of particles
N_list = [2:2:10;]
r = 10
d0_lb = 5e-10 # Minimum distance between the atoms (lower boundary) in m
λ = 421e-9
d0_lb = d0_lb / λ
# Quantization axis along z
e = [0, 0, 1.];

u_SS_N, p0_N = [], []

for (i, N) in enumerate(N_list)
    u_SS_r, p0_r = [], []
    for j in 1:r
        u_SS = stationnary_state(N, r, j)
        push!(u_SS_r, u_SS)

        @load "SS/r0/r0_N_$(N)_r_$j.jdl2" r0 L
        push!(p0_r, compute_params(N, r0, e))
    end
    push!(u_SS_N, u_SS_r), push!(p0_N, p0_r)
end

# # # Prepare the wrapper
# const N_FUNCS = length(N_list) # Total function nbr
# const functs = Vector{Function}(undef, N_FUNCS)

# for (i, N) in enumerate(N_list)
#     libpath ="Decay/libs/liballfuncs_$N.dll"
#     functs[i] = load_f("diffeqf", libpath)
# end

for (i, N) in enumerate(N_list)
    p0_r = []
    for j in 1:r
        if u_SS_N[i][j] != false
            du_SS = zeros(ComplexF64, length(u_SS_N[i][j]))
            libpath = r"Decay/libs/liballfuncs_$N.dll"
            @ccall(diffeqf, Cvoid, (Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}), du_SS, u_SS_N[i][j], p0_N[i][j])
            
            # functs[i](du_SS, u_SS_N[i][j], p0_N[i][j])
            println(du_SS)
        end
    end
end