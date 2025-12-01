# using IJulia
# installkernel("Julia (5 threads)", env=Dict("JULIA_NUM_THREADS"=>"5"))
# println(Threads.nthreads())

# To run multiple threads, we do not use local environment... -> Better solution to work on

import Pkg
Pkg.add(Pkg.PackageSpec(name="QuantumCumulants", version="0.3.8"))
Pkg.add("JLD2")
Pkg.add("PyPlot")
Pkg.add("Statistics")
Pkg.add("Symbolics")
Pkg.add("Unitful")
Pkg.add("PhysicalConstants")
Pkg.add("OrdinaryDiffEq")
Pkg.add("NonlinearSolve")
Pkg.add("ProgressBars")
Pkg.add("ProgressMeter")
Pkg.add("BenchmarkTools")