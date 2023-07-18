include("PinnBearing.jl")
using  Statistics, BenchmarkTools
using .PinnBearing 
using Plots

plotly()
nx = 380;ny = 241

bearing_lin = PinnBearing.HD_Bearing_Pfeil()
arc_vec = [64,4,128,4,64]
prob = PinnBearing.DNetPdeProblem(nx,ny,bearing_lin,arc_vec,"models/linear/model_parameters_DON6.hdf5")

state_vec = [0.98f0 * bearing_lin.c,0,0,0.0,0,0,pi/2,0]
alpha_wl = state_vec[5]; e = state_vec[1] / bearing_lin.c
eps = 2 * (sqrt(1 - e^2 * sin(alpha_wl)^2) -  e *abs(cos(alpha_wl))) * 0.95
state_vec[6] = eps
fx_lin,fy_lin,Mx,My = PinnBearing.bearing_pressure(state_vec,prob)

fx_dl,fy_dl,Mx_dl,My_dl = PinnBearing.forces_dl(prob,state_vec)
println("fx_fdm, fy_fdm: ",fx_lin," ",fy_lin," ", sqrt(fx_lin^2 + fy_lin^2))
println("fx_dl, fy_dl: ",fx_dl," ",fy_dl," ", sqrt(fx_dl^2 + fy_dl^2))
println("Mx_fdm, My_fdm: "," ",Mx,My)
println("Mx_dl, My_dl: ",Mx_dl," ",My_dl)
# arc_vec = [64,4,64,4,64]

# bearing = PinnBearing.Foil_Bearing_Leister()
# prob = PinnBearing.DNetPdeProblem_nl(nx,ny,bearing,arc_vec,"models/nonlinear/model_parameters_nl_DON3.hdf5")

# bearing.om = 400 * 2 * pi


# state_vec = [2.0f0 * bearing.c,0]
# fx,fy = PinnBearing.forces_dl_nl(prob,state_vec)

# println("fx: ",fx," fy: ",fy)

# fx2,fy2 = PinnBearing.nonlinear_pressure(state_vec,prob)

# println("fx2: ",fx2," fy2: ",fy2)







nothing