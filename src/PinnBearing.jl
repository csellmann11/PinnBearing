"""
PinnBearing.jl

Description:
This module contains the functions to calcutate the pressure and reaction forces
in a radial bearing. The forces can be calculated with a DeepONet or with a
finite difference method (FDM).
"""
module PinnBearing
    using LoopVectorization
    using MKL, Lux, Trapz, Plots, LoopVectorization, LinearAlgebra
    using KLU, BenchmarkTools
    include("DeepONet/DeepONet.jl")
    include("Physics/Bearing_comp/bearing.jl")
    include("Physics/PDE.jl")
    include("calc_forces.jl")
    include("calc_forces_nonlinear.jl")
    include("Physics/Bearing_comp/pressureFDM.jl")
    include("Physics/Bearing_comp/pressureSOR.jl")
    
    
    export HD_Bearing_Pfeil, HD_Bearing, HD_Bearing_Schweizer, Foil_Bearing_Leister, Foil_Bearing
    export DNetPdeProblem, PdeProb, forces_dl,DNetPdeProblem_nl, forces_dl_nl
    export DeepONet, FFCN
    export bearing_pressure
end 
