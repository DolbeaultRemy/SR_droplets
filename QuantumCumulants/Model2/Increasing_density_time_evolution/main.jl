import PhysicalConstants.CODATA2018: c_0
using QuantumCumulants
using CollectiveSpins
using Symbolics
using JLD2
using Unitful
using LaTeXStrings
using ProgressBars
using PyPlot

""" Create the big function that will call all the subfunctions to avoid doing a lot of ccalls """
function generate_full_dispatcher(filename::String, n::Int)
    open(filename, "w") do io
        println(io, "#include <complex.h>\n#include <math.h>\n")
        # External subfunctions
        for i in 1:(n)
            println(io, "extern void diffeqf_$i(double complex* du, const double complex* RHS1);")
        end

        println(io, "\nvoid diffeqf(double complex* du, const double complex* RHS1) {")

        # Call the subfunctions
        for i in 1:(n)
            println(io, "    diffeqf_$i(du, RHS1);")
        end

        println(io, "}")
    end
end



""" Create the object file (with the name of all the functions and their corresponding files) to avoid compilation issue """
function objs_file(eqs_eval)
    open("objs.txt", "w") do io
         println(io, "dispatcher.o")
        for i in 1:length(eqs_eval)
            println(io, "Cfunctions/diffeqf_$i.o")
        end
    end
end


""" Function creating the makefile for the corresponding HT frequency """
function change_mkfile(N)
    write("Makefile", replace(read("MakefileTemplate.txt", String), "liballfuncs.dll"=>"libs/liballfuncs_$N.dll"))
end

# Nbr of particles
N_list = [2:50;]

# Normalisation parameters
λ = 421e-9
γ = 32.7e6 # In Hz

# Physical values
ω0 = (2π*ustrip(c_0)/λ)
ωl = ω0
kl = [ustrip(c_0)/ωl, 0, 0] # Laser along x
Ω_RO = 1e7 # Taken from Barbut arXiv:2412.02541v1


# Position of atoms
Lx, Ly, Lz = [1, 1, 1] * 1e-6

# Normalization
ω0 = ω0 / γ
ωl = ωl / γ
kl = kl * λ
Ω_RO = Ω_RO / γ



# Quantization axis along z
e = [0, 0, 1.]

# Derivation of the symbolic MPC equations
@cnumbers Nsymb ΩROs

h = NLevelSpace(Symbol(:atom),2)

exp_RO(i) = IndexedVariable(:expRO, i)
conj_exp_RO(i) = IndexedVariable(:conjexpRO, i)
Γ(i,j) = IndexedVariable(:Γ,i,j)
Ω(i,j) = IndexedVariable(:Ω,i,j;identical=false)

i = Index(h,:i,Nsymb,h)
j = Index(h,:j,Nsymb,h)
k = Index(h,:k,Nsymb,h)

σ(x,y,z) = IndexedOperator(Transition(h,:σ,x,y),z)

H_RO = ΩROs * ∑(exp_RO(i)*σ(2,1,i) + conj_exp_RO(i)*σ(1,2,i), i)
H_elec = Σ(Ω(i,j)*σ(2,1,i)*σ(1,2,j), j, i)

H = Symbolics.simplify(H_RO + H_elec)

J = [σ(1,2,i)] # σ-, jump operators for the Lindbladian
rates = [Γ(i,j)]

ops = [σ(2, 2, k), σ(2, 1, k)]; # n_up/σ+

# Compute the MPC equations
eqs = meanfield(ops,H,J;rates=rates,order=2)
complete!(eqs);

# Create the directories
if !isdir("Image_positions")
    mkdir("Image_positions")
end
if !isdir("Cfunctions")
    mkdir("Cfunctions")
end
if !isdir("r0")
    mkdir("r0")
end
if !isdir("libs")
    mkdir("libs")
end
if !isdir("op_list")
    mkdir("op_list")
end

for N in ProgressBar(N_list)
    # Save the atoms position for comparison with QuantumOptics
    r0 = [[rand(Float64)*Lx, rand(Float64)*Ly, rand(Float64)*Lz] for i in 1:N]
    @save "r0/r0_N_$N.jdl2" r0
    r0 = r0 / λ

    # Plot the position of atoms
    plt.close("all")
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.scatter([r[1] for r in r0], [r[2] for r in r0], [r[3] for r in r0])
    ax.set_xlabel(L"x/$\lambda$")
    ax.set_ylabel(L"y/$\lambda$")
    ax.set_zlabel(L"z/$\lambda$")
    plt.savefig("Image_positions/Image_positions_N_$N")
    
    # Compute the Ω and Γ matrices of the electric dipole-dipole interaction using CollectiveSpins
    system = SpinCollection(r0, e, gammas=1.)
    Ω_CS = interaction.OmegaMatrix(system)
    Γ_CS = interaction.GammaMatrix(system)
    
    # Evaluate equations
    eqs_eval = evaluate(eqs; limits=(Nsymb=>N))

    # Evaluate parameters
    Γij_symb = [Γ(i,j) for i = 1:N for j=1:N]
    Ωij_symb = [Ω(i,j) for i = 1:N for j=1:N if i≠j]
    exp_RO_symb = [exp_RO(i) for i = 1:N]
    conj_exp_RO_symb = [conj_exp_RO(i) for i = 1:N]
    ps = [Γij_symb; Ωij_symb; exp_RO_symb; conj_exp_RO_symb; ΩROs]

    Γij_ = [Γ_CS[i, j] for i = 1:N for j=1:N]
    Ωij_ = [Ω_CS[i, j] for i = 1:N for j=1:N if i≠j]
    exp_RO_ = [exp(1im*r0[i]'kl) for i =1:N]
    conj_exp_RO_ = [exp(-1im*r0[i]'kl) for i =1:N]
    p0 = [Γij_; Ωij_; exp_RO_; conj_exp_RO_; Ω_RO/2]
    p0 = ps .=> p0
    eqs_eval = substitute(eqs_eval, Dict(p0))


    # Compute variables
    op_list = []
    var_array = []
    for i in 1:length(eqs_eval)
        var = eqs_eval[i].lhs
        push!(var_array, var)
        
        v_str = string(var)
        em = eachmatch(r"σ(\d+)", v_str)
        ind = [m.captures[1] for m in em]
        push!(op_list, [parse(Int, i) for i in ind])
    end
    
    @save "op_list/op_list_N_$N.jdl2" op_list

    # Build function
    for i in 1:length(eqs_eval)
        # Save the C function
        code = Symbolics.build_function([eqs_eval[i].rhs], var_array, target=Symbolics.CTarget())
        c_code = replace(code, "im" => "*I", "double* du" => "double complex* du", "const double* RHS1" => "const double complex* RHS1", "du[0]" => "du[$(i-1)]", "diffeqf" => "diffeqf_$i")
        open("Cfunctions/diffeqf_$i.c", "w") do io
            print(io, "#include <complex.h>\n") # Include complex package
            print(io, c_code)
        end
    end
    # Free RAM memomry
    code = nothing
    c_code = nothing

    # Build dispatcher
    generate_full_dispatcher("dispatcher.c", length(eqs_eval))

    # Build object file
    objs_file(eqs_eval)

    # Change the makefile
    change_mkfile(N)

    # Compile C functions
    run(`make -j5`)
end