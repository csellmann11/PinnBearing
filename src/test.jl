include("PinnBearing.jl")

using .PinnBearing 

nx = 150;ny = 60

bearing = PinnBearing.HD_Bearing_Schweizer()
arc_vec = [64,5,128,4,64]
prob = PinnBearing.DNetPdeProblem(nx,ny,bearing,arc_vec,"models/model_parameters_DON5.hdf5")


omega = 600 * 2 * pi; state_vec = [0.90 * bearing.c,0.0,0.0,0,0,0.0,0,0];
bearing.om = omega
fx1,fy1 = PinnBearing.forces_dl(prob,state_vec)

fx2,fy2 = PinnBearing.bearing_pressure(state_vec,prob)

println("fx1: ",fx1," fx2: ",fx2," fy1: ",fy1," fy2: ",fy2)