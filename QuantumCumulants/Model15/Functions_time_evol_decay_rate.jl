using PyPlot
using Statistics
using JLD2
using OrdinaryDiffEq
using DiffEqCallbacks
import PhysicalConstants.CODATA2018: c_0
using Unitful
using ProgressMeter
using ProgressBars
using NonlinearSolve
using BenchmarkTools
using Libdl
using LinearAlgebra
include("MiniCollectiveSpins.jl")


""" Check if the distribution worked, return the SS vector, else false """
function stationnary_state(N, r, i) # Nbr of atoms, nbr of repetitions, current repetition
    @load "SS/solutions/sol_N_$(N)_r_$r.jld2" sol_tasks
    nbr_error =  vcat([s[3] for s in sol_tasks]...)
    if (length(nbr_error)>0) && (i ∈ nbr_error)
        return false
    else
        sol_t = vcat([s[4] for s in sol_tasks]...)[i]
        return sol_t[end]
    end
end


""" Function loading the different libraries """
function load_f(fname::String, libpath::String)
	lib = Libdl.dlopen(libpath)
	fptr = Libdl.dlsym(lib, fname)
	return (du, u, params) -> ccall(fptr, Cvoid, (Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}), du, u, params)
end

""" Return the minimum distance of a distribution of atoms r0 """
function min_r0(r0)
    N = length(r0)
    d0 = zeros(N, N) # Repetiton, atom i, distance from atom j
    for j in 1:N
        for k = 1:N
            d0[j, k] = norm(r0[j]-r0[k])
        end
    end
    return minimum(d0[d0 .> 0])
end

function compute_params(N, r0, e)
    system = SpinCollection(r0, e, gammas=1.)
    Ω_CS = OmegaMatrix(system)
    Γ_CS = GammaMatrix(system)
    Γij_ = [Γ_CS[i, j] for i = 1:N for j=1:N]
    Ωij_ = [Ω_CS[i, j] for i = 1:N for j=1:N if i≠j]
    p0 = ComplexF64.([Γij_; Ωij_])
    return p0
end

# Create the directories
if !isdir("solutions")
    mkdir("solutions")
end