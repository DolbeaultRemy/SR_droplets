using QuantumCumulants
using Symbolics
using JLD2
using ProgressMeter
using ProgressBars

""" Create the big function that will call all the subfunctions to avoid doing a lot of ccalls """
function generate_full_dispatcher(filename::String, n::Int)
    open(filename, "w") do io
        println(io, "#include <complex.h>\n#include <math.h>\n")
        # External subfunctions
        for i in 1:(n)
            println(io, "extern void diffeqf_$i(double complex* du, const double complex* RHS1, const double complex* parameters);")
        end

        println(io, "\nvoid diffeqf(double complex* du, const double complex* RHS1, const double complex* parameters) {")

        # Call the subfunctions
        for i in 1:(n)
            println(io, "    diffeqf_$i(du, RHS1, parameters);")
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


""" Function that writes the Cfunctions, used to be distributed between all the threads """
function write_c_func(eqs, chunk, var_array, p0)
    for i in chunk
        code = Symbolics.build_function([eqs[i].rhs], var_array, target=Symbolics.CTarget())
        c_code = replace(code, 
            "im" => "*I", "double* du" => "double complex* du", "const double* RHS1" => "double complex* RHS1, const double complex* parameters", "du[0]" => "du[$(i-1)]", "diffeqf" => "diffeqf_$i")
        # Replace the parameters by a parametor vector
        c_code = replace(c_code, p0...)
        open("Cfunctions/diffeqf_$i.c", "w") do io
            print(io, "#include <complex.h>\n") # Include complex package
            print(io, c_code)
        end
    end
end

# Nbr of atoms
N = 2

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

eqs = meanfield(ops,H,J;rates=rates,order=2)
complete!(eqs);

# Create the directories
if !isdir("Cfunctions")
    mkdir("Cfunctions")
end
if !isdir("libs")
    mkdir("libs")
end
if !isdir("op_list")
    mkdir("op_list")
end

# Compute MPC equations
eqs_eval = QuantumCumulants.evaluate(eqs; limits=(Nsymb=>N))

# Replace parameters by a big vector
Γij_symb = [Γ(i,j) for i = 1:N for j=1:N]
Ωij_symb = [Ω(i,j) for i = 1:N for j=1:N if i≠j]
exp_RO_symb = [exp_RO(i) for i = N:-1:1] # We go in the decreasing direction to avoid exp_RO(10) being replace by exp_RO(1)0
conj_exp_RO_symb = [conj_exp_RO(i) for i = N:-1:1]
ps = [[string(Γ) for Γ in Γij_symb]; [string(Ω) for Ω in Ωij_symb]; [string(exp_RO) for exp_RO in exp_RO_symb]; [string(conj_exp_RO) for conj_exp_RO in conj_exp_RO_symb]; string(ΩROs)]

p0 = ["parameters[$(i-1)]" for i = 1:length(ps)]
p0 = ps .=> p0

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
@save "op_list/op_list_$N.jdl2" op_list

# Build Cfunctions in parallel
list_eqs = 1:length(eqs_eval)
chunks = Iterators.partition(list_eqs, cld(length(list_eqs), Threads.nthreads()))
tasks = map(chunks) do chunk # Split the equations in chunck solved on each core
    Threads.@spawn write_c_func(eqs_eval, chunk, var_array, p0)
end
wait.(tasks)

generate_full_dispatcher("dispatcher.c", length(eqs_eval))

# Build object file
objs_file(eqs_eval)

# Change the makefile
change_mkfile(N);