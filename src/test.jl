include("PinnBearing.jl")
using MKL
using .PinnBearing, Plots

include("Physics/pressureFDM.jl")

plotly()
nx = 60; ny = 30;

reduced_dim = 64
b_hidden, b_depth = 64, 4
branch_nets = [[1,reduced_dim,b_depth,b_hidden],[1,reduced_dim,b_depth,b_hidden]
        ,[2,reduced_dim,b_depth,b_hidden],[1,reduced_dim,b_depth,b_hidden],
        [2,reduced_dim,b_depth,b_hidden],[2,reduced_dim,b_depth,b_hidden],[1,reduced_dim,b_depth,b_hidden]]

t_hidden, t_depth = 128, 4
trunc_net = [3,reduced_dim,t_depth,t_hidden]

bearing = PinnBearing.HD_Bearing_Schweizer()
prob = PinnBearing.PDE_Prob(nx,ny,bearing,branch_nets,trunc_net,"models/model_parameters_misa_fullV2.hdf5")



Dm = 0.9                           
omega = 600 * 2 * pi; state_vec = [0.5 * bearing.c,0.0,1.0,0,0,0.0,20.0,10.0];
e = sqrt(state_vec[1]^2 + state_vec[2]^2)/bearing.c

alpha = state_vec[5]
eps_max = 2 * (sqrt(1 - e^2 * sin(alpha)^2) - e * abs(cos(alpha)))

eps_ = eps_max * Dm
state_vec[6] = eps_

bearing.om = omega
fx,fy = PinnBearing.forces_dl(prob,state_vec)

println("fx: ",fx," fy: ",fy)


(fxr,fyr),init_prob = bearing_pressure(state_vec,prob)

println("fx: ",fxr," fy: ",fyr)

rel_x = abs((fxr - fx)/fxr)*100; rel_y = abs((fyr - fy)/fyr)*100

println("rel_x: ",rel_x,"% rel_y: ",rel_y, "%")
nothing

