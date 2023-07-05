module PinnBearing
    using MKL, Lux, Trapz, Plots, LoopVectorization, LinearAlgebra

    include("DeepONet/DeepONet.jl")
    include("Physics/bearing.jl")
    include("Physics/PDE.jl")
    

    export HD_Bearing_Pfeil, HD_Bearing, HD_Bearing_Schweizer
    export PDE_Prob, forces_dl
    export DeepONet
end 
