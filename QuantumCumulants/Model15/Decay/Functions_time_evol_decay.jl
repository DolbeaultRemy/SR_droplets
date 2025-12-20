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
function stationnary_state(N, r, i)
    @load "../SS/solutions/sol_N_$(N)_r_$r.jld2" sol_tasks
    nbr_error =  vcat([s[3] for s in sol_tasks]...)
    # println("$N, $r, $i")
    if (length(nbr_error)>0) && (i ∈ nbr_error)
        return false
    else
        sol_t = vcat([s[4] for s in sol_tasks]...)[i]
        return sol_t[end]
    end
end


""" Load the corresponding distribution and computes the decay rate at t=0.
If compute_t_evolution, compute the whole evolution. """
function solve_decay_SS(chunk, functs, N, r, f_counter, save_sol)
    list_t, popup_t, nbr_error_t, sol_t = [], [], [], []

    for i in chunk
        @load "../SS/r0/r0_N_$(N)_r_$i.jdl2" r0 L

        system = SpinCollection(r0, e, gammas=1.)
        Ω_CS = OmegaMatrix(system)
        Γ_CS = GammaMatrix(system)
        Γij_ = [Γ_CS[i, j] for i = 1:N for j=1:N]
        Ωij_ = [Ω_CS[i, j] for i = 1:N for j=1:N if i≠j]
        p0 = ComplexF64.([Γij_; Ωij_])
        
        # Load the functions
        fsolve(du, u, p, t) = functs[f_counter](du, u, p0)

        u0 = stationnary_state(N, r, i)
        if u0 == false # If u0 did not CV
            continue
        end

        prob = OrdinaryDiffEq.ODEProblem(fsolve, u0, (0, t_end))
        
        sol = OrdinaryDiffEq.solve(prob, OrdinaryDiffEq.DP5(); abstol=1e-15, reltol=1e-15) # ; dtmin = 1e-2

        push!(list_t, sol.t)
        push!(popup_t, [sum(real(sol.u[i][1:N])) for i=1:length(sol.t)])

        if save_sol
            push!(sol_t, sol.u)
        end

        if ~SciMLBase.successful_retcode(sol) # | (sol.t[end] == t_end) | ((popup_t[end][end] < 0) | (popup_t[end][end] > 1)) # If solution does not CV, add it to the error
            push!(nbr_error_t, i)
        end
    end
    return list_t, popup_t, nbr_error_t, sol_t
end

""" Function loading the block subfunction when a lot of equations are involved """
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

# Create the directories
if !isdir("solutions")
    mkdir("solutions")
end