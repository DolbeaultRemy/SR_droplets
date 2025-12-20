include("MiniCollectiveSpins.jl")
include("Functions_time_evol_decay_rate.jl")
using PyPlot
using Statistics
using JLD2

### Define the system ###
# Nbr of particles
N_list = [2:2:10;]
r = 10
d0_lb = 5e-10 # Minimum distance between the atoms (lower boundary) in m
λ = 421e-9
d0_lb = d0_lb / λ

# Quantization axis along z
e = [0, 0, 1.];

# # Prepare the wrapper
# const N_FUNCS = length(N_list) # Total function nbr
# const functs = Vector{Function}(undef, N_FUNCS)

# for (i, N) in enumerate(N_list)
#     libpath ="Decay/libs/liballfuncs_$N.dll"
#     functs[i] = load_f("diffeqf", libpath)
# end

deriv_SS_N = []

N = 2
j = 1 
u_SS = stationnary_state(N, r, j)
du_SS = zeros(length(u_SS))
@load "SS/r0/r0_N_$(N)_r_$j.jdl2" r0 L
p0 = compute_params(N, r0, e)


for i = 1:5
    @ccall "Decay/libs/liballfuncs_2.dll".diffeqf(
            du_SS::Ptr{ComplexF64},
            u_SS::Ptr{ComplexF64},
            p0::Ptr{ComplexF64}
        )::Cvoid
        println(du_SS)
end
# for (i, N) in enumerate(N_list)
#     deriv_SS = []
#     for j in 1:r
#         u_SS = stationnary_state(N, r, j)

#         if u_SS != false
#             @load "SS/r0/r0_N_$(N)_r_$j.jdl2" r0 L
#             p0 = compute_params(N, r0, e)
#             du_SS = zeros(length(u_SS))

#             # lib = Libdl.dlopen("Decay/libs/liballfuncs_$N.dll") # Open the library explicitly.
#             # sym = Libdl.dlsym(lib, "diffeqf")   # Get a symbol for the function to call.
#             # @ccall sym(du_SS::Ptr{ComplexF64}, u_SS::Ptr{ComplexF64}, p0::Ptr{ComplexF64})::Cvoid#, Cvoid, (Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}), )
#             # Libdl.dlclose(lib) # Close the library explicitly.
#             @ccall "Decay/libs/liballfuncs_2.dll".diffeqf(
#                     du_SS::Ptr{ComplexF64},
#                     u_SS::Ptr{ComplexF64},
#                     p0::Ptr{ComplexF64}
#                 )::Cvoid
#             # functs[i](du_SS, u_SS, p0)
#             # f(du_SS, u_SS)
#             # println(du_SS)
#             # # Libdl.dlclose(lib) # Close the library explicitly.

#             # push!(deriv_SS, du_SS)
#             # println(@ccall clock()::Int32)
#             # @ccall function_name(argvalue1::argtype1, ...)::returntype
#         end
#     end
#     # @save "solutions/derivative_N_$(N)_r_$(r).jld2" deriv_SS_N
# end
print("done")