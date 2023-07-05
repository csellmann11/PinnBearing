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

