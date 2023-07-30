module PinnBearing
    using MKL, Lux, Trapz, Plots, LoopVectorization, LinearAlgebra
    using KLU, BenchmarkTools
    include("DeepONet/DeepONet.jl")
    include("Physics/Bearing_comp/bearing.jl")
    include("Physics/PDE.jl")
    include("calc_forces.jl")
    include("calc_forces_nonlinear.jl")
    include("Physics/Bearing_comp/pressureFDM.jl")
    include("Physics/Bearing_comp/pressureSOR.jl")

    export HD_Bearing_Pfeil, HD_Bearing, HD_Bearing_Schweizer, Foil_Bearing_Leister
    export DNetPdeProblem, PdeProb, forces_dl
    export DeepONet, FFCN
    export bearing_pressure
end 
