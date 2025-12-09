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


""" Prepare the initial vector u0 """
function u0_CFunction(phi_array, theta_array, op_list)
    u0 = ones(ComplexF64, length(op_list))
    for i in 1:length(op_list)
        if length(op_list[i]) == 1
            j = Int(op_list[i][1] % 10^floor(log10(abs(op_list[i][1]))-1)) # Atom nbr
            if parse(Int, string(op_list[i][1])[1:2]) == 22
                u0[i] = cos(theta_array[j]/2)^2
            elseif parse(Int, string(op_list[i][1])[1:2]) == 21
                u0[i] = cos(theta_array[j]/2)*exp(1im*phi_array[j])*sin(theta_array[j]/2)
            else
                println(op_list[i][1])
            end
        end

        if length(op_list[i]) == 2
            for op in op_list[i]
                j = Int(op % 10^floor(log10(abs(op))-1)) # Atom nbr
                if parse(Int, string(op)[1:2]) == 22
                    u0[i] *= cos(theta_array[j]/2)^2
                elseif parse(Int, string(op)[1:2]) == 21
                    u0[i] *= cos(theta_array[j]/2)*exp(1im*phi_array[j])*sin(theta_array[j]/2)
                elseif parse(Int, string(op)[1:2]) == 12
                    u0[i] *= cos(theta_array[j]/2)*exp(-1im*phi_array[j])*sin(theta_array[j]/2)
                else
                    println(op)
                end
            end
        end
    end
    return u0
end


""" Create a random distribution, save it, computes the corresponding parameters an return the stationnary state. 
If compute_t_evolution, compute the whole evolution, else only the stationnary state. """
function solve_random_distrib(chunk, f, op_list, N, n, d0_lb, window_t, window_var, threshold_box, f_counter, save_sol)
    list_t, popup_t, nbr_error_t, sol_t = [], [], [], []

    for i in chunk
        # Compute distribution
        L = (N/n)^(1/3) # Change the volume to keep the density cste
        r0 = [[rand(Float64)*L, rand(Float64)*L, rand(Float64)*L] for i in 1:N]

        # Choose a distribution where the minimum distance between the atoms is bigger than d0_min
        while min_r0(r0) < d0_lb
            r0 = [[rand(Float64)*L, rand(Float64)*L, rand(Float64)*L] for i in 1:N]
        end

        # Save the atoms position for comparison with QuantumOptics
        @save "r0/r0_N_$(N)_r_$i.jdl2" r0 L

        # Compute the parameters
        system = SpinCollection(r0, e, gammas=1.)
        Ω_CS = OmegaMatrix(system)
        Γ_CS = GammaMatrix(system)
        Γij_ = [Γ_CS[i, j] for i = 1:N for j=1:N]
        Ωij_ = [Ω_CS[i, j] for i = 1:N for j=1:N if i≠j]
        exp_RO_ = [exp(1im*r0[i]'kl) for i = N:-1:1] # We go in the decreasing direction to avoid exp_RO(10) being replace by exp_RO(1)0
        conj_exp_RO_ = [exp(-1im*r0[i]'kl) for i = N:-1:1]
        p0 = ComplexF64.([Γij_; Ωij_; exp_RO_; conj_exp_RO_; Ω_RO/2])
        
        # Load the functions
        fsolve(du, u, p, t) = functs[f_counter](du, u, p0)

        phi_array_0, theta_array_0 = zeros(N), ones(N)*π # We start from all the atoms in the GS
        u0 = u0_CFunction(phi_array_0, theta_array_0, op_list)

        function stop_condition(u, t, integrator)
            T = integrator.sol.t
            if (t > 2*window_var) && (length(T) > 1) && (t % window_var < t-T[end])                
                popup = [mean(real([integrator.sol.u[j][i] for i = 1:N])) for j = 1:length(T)]
                popup_smooth = runmean_window_t(popup, T, window_t)
                boxes = boxmean_window_t(popup_smooth, T[T .> window_var], window_var)
                if length(boxes) > 1
                    Delta_boxes, T_boxes = boxes[2:end] - boxes[1:end-1], [i*window_var for i in 1:div(T[end], window_var)-1]
                    if abs(Delta_boxes[end]) < threshold_box
                        return true
                    end
                end
                # Sanity check
                if (popup[end] < 0) | (popup[end] > 1)
                    return true
                end
            end
            return false
        end

        function affect!(integrator)
            terminate!(integrator) # If stop conditions are met, stop the solver
        end

        stop_cb = DiscreteCallback(stop_condition, affect!)

        prob = OrdinaryDiffEq.ODEProblem(fsolve, u0, (0, t_end))
        
        sol = OrdinaryDiffEq.solve(prob, OrdinaryDiffEq.DP5();
                    reltol=1.0e-5,
                    abstol=1.0e-5,
                    callback = stop_cb,
                    dtmin = 1e-2)

        push!(list_t, sol.t)
        push!(popup_t, [sum(real(sol.u[i][1:N])) for i=1:length(sol.t)])

        if save_sol
            push!(sol_t, sol.u)
        end

        if ~SciMLBase.successful_retcode(sol) | (sol.t[end] == t_end) # If solution does not CV, add it to the error
            push!(nbr_error_t, i)
        end
    end
    return list_t, popup_t, nbr_error_t, sol_t
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

""" Running average, with spacing one, and no zero padding for the first values (they are skipped) """
function runmean_window_t(array, T, window_t)
    mean_array = []
    for i = 1:length(array)
        if T[i] > window_t
            push!(mean_array, mean(array[abs.(T .- T[i]) .< window_t]))
        end
    end
    return mean_array
end

""" Create boxes that do not overlap, return the mean on each one of these boxes """
function boxmean_window_t(array, T, window_t)
    boxes_mean = []
    i = 1
    while true
        t_start_window = T[i]
        array_box = []
        while T[i] - t_start_window < window_t
            push!(array_box, array[i])
            if i < length(T)
                i += 1
            else
                push!(boxes_mean, mean(array_box))
                return boxes_mean
            end
        end
        push!(boxes_mean, mean(array_box))
    end
end


""" Function loading the block subfunction when a lot of equations are involved """
function load_f(fname::String, libpath::String)
	lib = Libdl.dlopen(libpath)
	fptr = Libdl.dlsym(lib, fname)
	return (du, u, params) -> ccall(fptr, Cvoid, (Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}), du, u, params)
end

# Create the directories
if !isdir("r0")
    mkdir("r0")
end
if !isdir("Images_distribution")
    mkdir("Images_distribution")
end
if !isdir("solutions")
    mkdir("solutions")
end