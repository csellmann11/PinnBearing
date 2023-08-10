include("../src/PinnBearing.jl")
using .PinnBearing

# number of grid points in x and y direction
nx = 80;ny = 30

# load example bearing from file
bearing = PinnBearing.HD_Bearing_Pfeil()
# file of trained model
filename = "models/linear/model_parameters_DON2.hdf5"
# arc vector for the DeepONet
arc_vec = [64,4,128,4,64]

# create PDE problem
pde_prob = PinnBearing.DNetPdeProblem(nx,ny,bearing,arc_vec,filename)

# create state vector for parallel clearance
state_vector = [bearing.c * 0.2,0,0,0] .|> Float32

fun = (u) -> PinnBearing.forces_dl(u,pde_prob)

fun2 = (u) -> PinnBearing.bearing_pressure(u,pde_prob)

# forces with PINN
fx,fy,Mx,My = fun(state_vector) 

# forces with FDM
fx2,fy2,Mx2,My2 = fun2(state_vector) 

nothing