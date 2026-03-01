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

""" Create a random distribution, save it, computes the corresponding parameters """
function create_random_distribution(N, n, d0_lb, Ω_RO)
    # Compute distribution
    L = (N/n)^(1/2) # Change the volume to keep the 2D density = n and the atom nbr = N
    r0 = [[rand(Float64)*L, rand(Float64)*L, rand(Float64)*L] for atm in 1:N]

    # Choose a distribution where the minimum distance between the atoms is bigger than d0_min
    # println("$(min_r0(r0)), $(d0_lb)")
    while min_r0(r0) < d0_lb
        r0 = [[rand(Float64)*L, rand(Float64)*L, rand(Float64)*L] for a in 1:N]
    end

    # Save the atoms position for comparison with QuantumOptics
    @save "r0/r0_N_$(N)_n_$n.jdl2" r0 L

    # Computes the parameter of this system
    system = SpinCollection(r0, e, gammas=1.)
    Ω_CS = OmegaMatrix(system)
    Γ_CS = GammaMatrix(system)
    Γij_ = [Γ_CS[i, j] for i = 1:N for j=1:N]
    Ωij_ = [Ω_CS[i, j] for i = 1:N for j=1:N if i≠j]
    exp_RO_ = [exp(1im*r0[i]'kl) for i = N:-1:1] # We go in the decreasing direction to avoid exp_RO(10) being replace by exp_RO(1)0
    conj_exp_RO_ = [exp(-1im*r0[i]'kl) for i = N:-1:1]
    return ComplexF64.([Γij_; Ωij_; exp_RO_; conj_exp_RO_; Ω_RO])
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

""" For one chunck of density, computes the parameters, solves the system until it converges, returns [list_t, popup_t, sol_t].
If compute_t_evolution, returns the whole evolution, else only the SS. """
function solve_random_distrib_vary_Ω(chunk, fsolve_diff_eq, op_list, N, n_list, d0_lb, window_t, window_var, threshold_box, t_end, save_sol, Ω_RO)
    list_t, popup_t, sol_t = [], [], []

    # Initial state
    phi_array_0, theta_array_0 = zeros(N), ones(N)*π # We start from all the atoms in the GS
    u0 = u0_CFunction(phi_array_0, theta_array_0, op_list)

    # Defining the stoping functions (in solve fct to handle the global variables [window_t, window_var, threshold_box])
    """ Function that tells the solver when it is possible to stop (no more evolution of the system) """
    function stop_condition(u, t, integrator)
        T = integrator.sol.t
        if (t > 2*window_var) && (length(T) > 1) && (t % window_var < t-T[end]) && (length(T[T .> window_var]) > 0)               
            popup = [mean(real([integrator.sol.u[j][i] for i = 1:N])) for j = 1:length(T)]
            popup_smooth = runmean_window_t(popup, T, window_t)
            boxes = boxmean_window_t(popup_smooth, T[T .> window_var], window_var)
            if length(boxes) > 1
                Delta_boxes = boxes[2:end] - boxes[1:end-1]
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

    """ Function to terminate the solver when stop_condition() is fulfilled """
    function affect!(integrator)
        terminate!(integrator) # If stop conditions are met, stop the solver
    end

    stop_cb = DiscreteCallback(stop_condition, affect!)

    for i in chunk
        CV_sim = false # Criterium checking if simulations converged
        sol = nothing # Define sol for multithreading

        try
            while !CV_sim
                p0 = create_random_distribution(N, n_list[i], d0_lb, Ω_RO) # Parameters of the random distribution of density n_list[i]
                
                fsolve(du, u, p, t) = fsolve_diff_eq(du, u, p0)

                # Define the ODE
                prob = OrdinaryDiffEq.ODEProblem(fsolve, u0, (0, t_end))
                sol = OrdinaryDiffEq.solve(prob, OrdinaryDiffEq.DP5(); callback = stop_cb)
                
                Nup_t = [sum(real(sol.u[i][1:N])) for i=1:length(sol.t)] # Expectation values of the number of atoms in the ES

                if SciMLBase.successful_retcode(sol) && (sol.t[end] < t_end) && (minimum(Nup_t) > 0) && (maximum(Nup_t) < N) # If solution CV and is physical, save the solution
                    CV_sim = true
                else
                    println("For a normalised density of $(n_list[i]), the solver couldn't CV... Lauching a new simulation.")
                    println("tend = $(sol.t[end])")
                    println("max(Nup) = $(maximum(Nup_t))")
                    println("min(Nup) = $(minimum(Nup_t))")
                    println()
                end
            end

        catch err
            println("Simulation failed for normalised density = $(n_list[i])")
            println(err)
        end

        @assert sol !== nothing # Sanity check to see if the simulation was actually performed

        # Save solutions
        push!(list_t, sol.t)
        push!(popup_t, [sum(real(sol.u[i][1:N])) for i=1:length(sol.t)])

        if save_sol # Save the entire set of expectation values w.r.t. time
            push!(sol_t, sol.u)
        else
            push!(sol_t, sol.u[end]) # Only keep the last state of the system if save_sol = false
        end
    end
    
    return list_t, popup_t, sol_t
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

# Create the directories
if !isdir("r0")
    mkdir("r0")
end

if !isdir("solutions")
    mkdir("solutions")
end