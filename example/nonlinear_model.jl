using BenchmarkTools
include("../src/PinnBearing.jl")
using .PinnBearing

# number of grid points in x and y direction
nx = 60;ny = 30

# load example bearing from file
bearing = PinnBearing.Foil_Bearing_Leister()
bearing.om = 400 * 2pi
# file of trained model
filename = "models/nonlinear/model_parameters_nl_DON2.hdf5"
# arc vector for the DeepONet
arc_vec = [32,4,64,4,32]

# create PDE problem
pde_prob = PinnBearing.DNetPdeProblem(nx,ny,bearing,arc_vec,filename)

# create state vector for parallel clearance and static state
state_vector = [bearing.c * 1.5,0] .|> Float32

fun = (u) -> PinnBearing.forces_dl_nl(u,pde_prob)

fun2 = (u) -> PinnBearing.nonlinear_bearing_pressure(u,pde_prob)

# forces with PINN
fx,fy,Mx,My = fun(state_vector) 

# forces with FDM
fx2,fy2,Mx2,My2 = fun2(state_vector)
nothing