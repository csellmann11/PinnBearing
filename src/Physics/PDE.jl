include("load_file.jl")
struct PDE_Prob
    nx :: Int
    ny :: Int
    dx :: Float32
    dy :: Float32
    
    x :: Array{Float32}
    y :: Array{Float32}
    X :: Array{Float32}
    Y :: Array{Float32}
    H :: Array{Float32}
    
    cosX    :: Array{Float32}
    sinX    :: Array{Float32}
    pressure  :: Array{Float32}

    bearing  :: HD_Bearing
    

    model       ::DeepONet
    model_pars  ::NamedTuple
    model_state ::NamedTuple
end

function PDE_Prob(nx::Int,ny::Int,bearing::HD_Bearing,branch_networks_info,trunc_network_info,
    params_file)
    dx :: Float32 = 2pi/(nx-1)
    dy :: Float32 = bearing.B/(ny-1)

    x = collect(LinRange(0,2pi,nx)) .|> Float32
    y = collect(LinRange(-1,1,ny)) .|> Float32
    
    X = [xv for xv in x, _ in y]
    Y = [yv for _ in x, yv in y]


    H, pressure = similar(X), similar(X)
    cosX = cos.(X); sinX = sin.(X)


    y_data = reshape(Y,1,nx*ny); x_data = reshape(X,1,nx*ny)

    model = DeepONet(branch_networks_info,trunc_network_info,y_data)
    model.flags.training = false
    rng = MersenneTwister(1234)

    ps,st = Lux.setup(rng,model)

    trunc_in = cat(cos.(x_data),sin.(x_data),y_data,dims=1) .|> Float32
    branch_in = [[0.0f0],[0.0f0],[0.0f0,0.0f0],[0.0f0],[0.0f0,0.0f0],[0.0f0,0.0f0],[0.0f0]]
    input = [trunc_in,branch_in...]

    OpenHDF5(params_file,ps)
    model(input,ps,st)

    PDE_Prob(nx,ny,dx,dy,x,y,X,Y,H,cosX,sinX,pressure,bearing,model,ps,st)
end

function forces_dl(prob::PDE_Prob,state_vector,ω)
    model, ps, st = prob.model, prob.model_pars, prob. model_state
    bearing = prob.bearing

    
    xs,ys,dxs,dys,alpha_WL,eps_,eps_dot,alpha_WL_dot = state_vector .|> Float32
    ω = Float32(ω); u_m = bearing.rI * ω/2

    E ::Float32 = sqrt(xs^2 + ys^2)/bearing.c
    Φ ::Float32 = atan(ys,xs) 

    _ϕ ::Float32 = (dys*xs - dxs*ys) / (xs^2 + ys^2)
    _e ::Float32 = (dxs*xs + dys*ys) / sqrt(xs^2 + ys^2)

    _Φ ::Float32 = _ϕ * bearing.rI/u_m
    _E ::Float32 = _e * bearing.rI/u_m/bearing.c 

    eps_max :: Float32 = 2 * (sqrt(1 - E^2*sin(alpha_WL)^2) - E * abs(cos(alpha_WL)))
    D_m = eps_/eps_max

    ϕ01 = atan((E * (_Φ-1)),_E)
    ϕ23 = atan(-D_m,0) ##TODO: Anapssen wenn Geschwindigkeit nicht konstant 0

    alpha_WL_in = [cos(alpha_WL),sin(alpha_WL)]
    ϕ01_in :: Vector{Float32} = [cos(ϕ01),sin(ϕ01)]
    ϕ23_in :: Vector{Float32} = [cos(ϕ23),sin(ϕ23)]

    const1 = -D_m * eps_max

    λ::Float32 = sqrt(_E^2 + (E * (_Φ-1))^2) + sqrt(const1^2)
    
    in1,in2,in3,in4,in5,in6 = [2*E-1], [2*sqrt(_E^2 + (E * (_Φ-1))^2)/λ - 1], ϕ01_in,[2*D_m - 1],alpha_WL_in,ϕ23_in
    width = 2; width_in = [2/3 * (width - 2.5)] .|> Float32

    inputs::Vector{Union{Vector{Float32},Float32}}   = [in1,in1,in2,in3,in4,in5,in6, width_in]

    p_fak = λ* bearing.eta * u_m * bearing.rI^2/
        (bearing.c^2) * bearing.b * pi/2


    @. prob.H = 1 + E*prob.cosX + D_m * eps_max * 1/2 * prob.Y * cos(prob.X - alpha_WL)
    model(inputs,ps,st)
    prob.pressure .= reshape(model.output,prob.nx,prob.ny)

    p_fak = λ* bearing.eta * u_m * bearing.rI^2/
        (bearing.c^2) * bearing.b/2
    @. prob.pressure = max(prob.pressure/prob.H^2,0) 

    fx = trapz((prob.x,prob.y),prob.pressure .* prob.cosX) * p_fak
    fy = trapz((prob.x,prob.y),prob.pressure .* prob.sinX) * p_fak

    Rot = [cos(Φ) -sin(Φ); sin(Φ) cos(Φ)]
    F = Rot * [fx;fy]
    fx,fy = F[1],F[2]

    return fx,fy
end