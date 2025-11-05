import Pkg
Pkg.activate(".")
Pkg.add(Pkg.PackageSpec(name="QuantumCumulants", version="0.3.8"))
Pkg.add(Pkg.PackageSpec(name="OrdinaryDiffEq", version="6.98.0"))
Pkg.add(Pkg.PackageSpec(name="Symbolics", version="6.31.1"))
Pkg.add(Pkg.PackageSpec(name="QuantumOptics", version="1.2.3")) # Has to be 1.2.4 to use timeevolution.master! Conflict with CollectiveSpins.
Pkg.add("PhysicalConstants")
Pkg.add("Unitful")
Pkg.add("PyPlot")
#Pkg.develop("CollectiveSpins")