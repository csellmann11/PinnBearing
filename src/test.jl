include("PinnBearing.jl")
using  Statistics, BenchmarkTools
using .PinnBearing 
using Plots

plotly()
# nx = 380;ny = 241

# bearing_lin = PinnBearing.HD_Bearing_Pfeil()
# arc_vec = [64,4,128,4,64]
# prob = PinnBearing.DNetPdeProblem(nx,ny,bearing_lin,arc_vec,"models/linear/model_parameters_DON6.hdf5")

# state_vec = [0.98f0 * bearing_lin.c,0,0,0.0,0,0,pi/2,0]
# alpha_wl = state_vec[5]; e = state_vec[1] / bearing_lin.c
# eps = 2 * (sqrt(1 - e^2 * sin(alpha_wl)^2) -  e *abs(cos(alpha_wl))) * 0.95
# state_vec[6] = eps
# fx_lin,fy_lin,Mx,My = PinnBearing.bearing_pressure(state_vec,prob)

# fx_dl,fy_dl,Mx_dl,My_dl = PinnBearing.forces_dl(prob,state_vec)
# println("fx_fdm, fy_fdm: ",fx_lin," ",fy_lin," ", sqrt(fx_lin^2 + fy_lin^2))
# println("fx_dl, fy_dl: ",fx_dl," ",fy_dl," ", sqrt(fx_dl^2 + fy_dl^2))
# println("Mx_fdm, My_fdm: "," ",Mx,My)
# println("Mx_dl, My_dl: ",Mx_dl," ",My_dl)

nx = 60; ny = 30
arc_vec = [32,4,64,4,32]

grid = [[25,10],[40,15],[60,25],[80,30],[100,40],[120,50]]#,[140,60]]
#grid = [[80,25],[80,35],[80,45],[80,55],[80,65]]
#grid = [[50,25],[70,25],[90,25],[110,25],[130,25]]
forces = Array{Float32}(undef,length(grid),2)   

bearing = PinnBearing.Foil_Bearing_Leister()

for i in eachindex(grid)
    prob = PinnBearing.DNetPdeProblem_nl(grid[i][1],grid[i][2],bearing,arc_vec,"models/nonlinear/model_parameters_nl_DON2.hdf5")

    bearing.om = 400 * 2 * pi


    state_vec = [1.9f0 * bearing.c,0]
    fx,fy = PinnBearing.forces_dl_nl(prob,state_vec)
    f = sqrt(fx^2 + fy^2)

    println("fx: ",fx," fy: ",fy)

    fx2,fy2 = PinnBearing.nonlinear_pressure(state_vec,prob)
    f2 = sqrt(fx2^2 + fy2^2)
    println("fx2: ",fx2," fy2: ",fy2)

    forces[i,:] = [f,f2]
end


x = 1:size(forces,1)
dl_last = forces[end,1]; fdm_last = forces[end,2]
dl_upper = 1.01 * dl_last * ones(size(forces,1)); dl_lower = 0.99 * dl_last * ones(size(forces,1))
fdm_upper = 1.01 * fdm_last * ones(size(forces,1)); fdm_lower = 0.99 * fdm_last * ones(size(forces,1))

p = plot(x,forces[:,1],label="PINN", linewidth = 2, color = :red)

plot!(x, dl_last * ones(size(forces,1)), fill = (dl_lower, 0.2, :red), linealpha = 0.0, color = :red, label = false)
plot!(x, dl_last * ones(size(forces,1)), fill = (dl_upper, 0.2, :red), linealpha = 0.0, color = :red, label = false)


plot!(x,forces[:,2],label="FDM", linewidth = 2, color = :green)

plot!(x, fdm_last * ones(size(forces,1)), fill = (fdm_lower, 0.2, :green), linealpha = 0.0, color = :green, label = false)
plot!(x, fdm_last * ones(size(forces,1)), fill = (fdm_upper, 0.2, :green), linealpha = 0.0, color = :green, label = false)


gr()
xticks!(1:length(grid),["25 x 10","40 x 15","60 x 25","80 x 30","100 x 40","120 x 50"])
xlabel!("Diskretisierung"), ylabel!("||F||[N]")
plot!(tickfont = font(10),
xlabelfontsize = 15,
ylabelfontsize = 15)
display(plot!())

savefig(p,"convergeny_SOR.pdf")


# plot!(x, dl_upper, label="", color=:red)
# plot!(x, dl_lower, color=:red)
# plot!(x, dl_upper, fill=(x, dl_lower), fillalpha=0.2, fillcolor=:red, label="DL ±1%")









nothing