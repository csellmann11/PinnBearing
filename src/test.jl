include("PinnBearing.jl")

using .PinnBearing
using BenchmarkTools


nx = 80;ny = 30;

arc_vec = [64,4,128,4,64]
filename = "models/linear/model_parameters_DON2.hdf5"

bearing = PinnBearing.HD_Bearing_Pfeil()
bearing.om = 200.0f0
prob = PinnBearing.DNetPdeProblem(nx,ny,bearing,arc_vec,filename)

state_vec = [0.9 * bearing.c,0,0.0,0] .|> Float32

@time for i in 1:1000
    PinnBearing.forces_dl(prob,state_vec, parallel = true)
end

PinnBearing.bearing_pressure(state_vec,prob,parallel = true)
@time for i in 1:1000
    PinnBearing.bearing_pressure(state_vec,prob,parallel = true)
end