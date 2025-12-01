import Pkg
Pkg.activate(".") # Change with working directory if you are not already inside
Pkg.add(Pkg.PackageSpec(name="QuantumCumulants", version="0.3.8"))
Pkg.add("JLD2")
Pkg.add("PyPlot")
Pkg.add("Statistics")
Pkg.add("Symbolics")
Pkg.add("Unitful")
Pkg.add("PhysicalConstants")
Pkg.add("SteadyStateDiffEq")
Pkg.add("NonlinearSolve")
Pkg.add("ProgressBars")
Pkg.add("ProgressMeter")

# CollectiveSpins is outdated. Import it in dev mode, then go in you julia/dev directory and change the package by the one of this github (after unzipping it)
#Pkg.develop("CollectiveSpins")
#Pkg.add(Pkg.PackageSpec(name="SteadyStateDiffEq", version="2.8.0"))