include("PinnBearing.jl")
using MKL
using .PinnBearing, Plots

plotly()
nx = 80; ny = 50;

reduced_dim = 64
b_hidden, b_depth = 64, 4
branch_nets = [[1,reduced_dim,b_depth,b_hidden],[1,reduced_dim,b_depth,b_hidden]
        ,[2,reduced_dim,b_depth,b_hidden],[1,reduced_dim,b_depth,b_hidden],
        [2,reduced_dim,b_depth,b_hidden],[2,reduced_dim,b_depth,b_hidden],[1,reduced_dim,b_depth,b_hidden]]

t_hidden, t_depth = 128, 4
trunc_net = [3,reduced_dim,t_depth,t_hidden]

bearing = PinnBearing.HD_Bearing_Pfeil()
prob = PinnBearing.PDE_Prob(nx,ny,bearing,branch_nets,trunc_net,"models/model_parameters_misa_fullV2.hdf5")



Dm = 0.0                           
omega = 600 * 2 * pi; state_vec = [0.9 * bearing.c,0.0,0.0,0,0,0.0,0,0];
e = sqrt(state_vec[1]^2 + state_vec[2]^2)/bearing.c

alpha = state_vec[5]
eps_max = 2 * (sqrt(1 - e^2 * sin(alpha)^2) - e * abs(cos(alpha)))

eps_ = eps_max * Dm
state_vec[6] = eps_

@time for i in 1:1000
        fx,fy = PinnBearing.forces_dl(prob,state_vec,omega)
end
println("fx: ",fx," fy: ",fy)


