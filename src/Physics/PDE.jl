
abstract type AbstractPdeProblem end 
abstract type AbstractDNetPdeProblem <: AbstractPdeProblem end


include("load_file.jl")
struct DNetPdeProblem <: AbstractDNetPdeProblem
    nx :: Int
    ny :: Int
    dx :: Float32
    dy :: Float32
    
    x :: Array{Float32}
    y :: Array{Float32}
    X :: Array{Float32}
    Y :: Array{Float32}
    #H :: Array{Float32}
    
    cosX    :: Array{Float32}
    sinX    :: Array{Float32}
    # pressure  :: Array{Float32}

    bearing  :: Bearing
    

    model       ::DeepONet
    model_pars  ::NamedTuple
    model_state ::NamedTuple
end

function DNetPdeProblem(nx::Int,ny::Int,bearing::HD_Bearing,arc_vec,
    params_file)
    dx :: Float32 = 2pi/(nx-1)
    dy :: Float32 = 1/(ny-1)

    x = collect(LinRange(0,2pi,nx)) .|> Float32
    y = collect(LinRange(-1,1,ny)) .|> Float32
    
    X = [xv for xv in x, _ in y]
    Y = [yv for _ in x, yv in y]


    #H, pressure = similar(X), similar(X)
    cosX = cos.(X); sinX = sin.(X)


    y_data = reshape(Y,1,nx*ny); x_data = reshape(X,1,nx*ny)

    b_hidden, b_depth = arc_vec[1], arc_vec[2]
    t_hidden, t_depth = arc_vec[3], arc_vec[4]
    reduced_dim = arc_vec[5]

    println("b_hidden: ",b_hidden," b_depth: ",b_depth," t_hidden: ",t_hidden," t_depth: ",t_depth," reduced_dim: ",reduced_dim)

    branch_networks_info = [[1,reduced_dim,b_depth,b_hidden],[1,reduced_dim,b_depth,b_hidden]
            ,[2,reduced_dim,b_depth,b_hidden],[1,reduced_dim,b_depth,b_hidden],
            [2,reduced_dim,b_depth,b_hidden],[2,reduced_dim,b_depth,b_hidden],[1,reduced_dim,b_depth,b_hidden]]

    trunc_network_info = [3,reduced_dim,t_depth,t_hidden]

    model = DeepONet(branch_networks_info,trunc_network_info,y_data)
    model.flags.training = false
    rng = MersenneTwister(1234)

    ps,st = Lux.setup(rng,model)

    trunc_in = cat(cos.(x_data),sin.(x_data),y_data,dims=1) .|> Float32
    branch_in = [[0.0f0],[0.0f0],[0.0f0,0.0f0],[0.0f0],[0.0f0,0.0f0],[0.0f0,0.0f0],[0.0f0]]
    input = [trunc_in,branch_in...]

    OpenHDF5(params_file,ps)
    model(input,ps,st)

    DNetPdeProblem(nx,ny,dx,dy,x,y,X,Y,cosX,sinX,bearing,model,ps,st)
end

function DNetPdeProblem_nl(nx::Int,ny::Int,bearing::Foil_Bearing,arc_vec,
    params_file; misa = false)
    dx :: Float32 = 2pi/(nx-1)
    dy :: Float32 = 1/(ny-1)

    x = collect(LinRange(0,2pi,nx)) .|> Float32
    y = collect(LinRange(-1,1,ny)) .|> Float32
    
    X = [xv for xv in x, _ in y]
    Y = [yv for _ in x, yv in y]


    #H, pressure = similar(X), similar(X)
    cosX = cos.(X); sinX = sin.(X)


    y_data = reshape(Y,1,nx*ny); x_data = reshape(X,1,nx*ny)

    b_hidden, b_depth = arc_vec[1], arc_vec[2]
    t_hidden, t_depth = arc_vec[3], arc_vec[4]
    reduced_dim = arc_vec[5]

    println("b_hidden: ",b_hidden," b_depth: ",b_depth," t_hidden: ",t_hidden," t_depth: ",t_depth," reduced_dim: ",reduced_dim)

    if misa
        branch_networks_info = [[1,reduced_dim,b_depth,b_hidden],[1,reduced_dim,b_depth,b_hidden],[1,reduced_dim,b_depth,b_hidden]
                ,[2,reduced_dim,b_depth,b_hidden]]

        branch_in = [[0.0f0],[0.0f0],[0.0f0],[0.0f0,0.0f0]]

        trunc_in = cat(cos.(x_data),sin.(x_data),y_data,dims=1) .|> Float32
    else
        branch_networks_info = [[1,reduced_dim,b_depth,b_hidden],[1,reduced_dim,b_depth,b_hidden]]
        branch_in = [[0.0f0],[0.0f0]]

        trunc_in = cat(cos.(x_data),sin.(x_data),cos.(y_data/2 * pi),dims=1) .|> Float32
    end
    trunc_network_info = [3,reduced_dim,t_depth,t_hidden]

    model = DeepONet(branch_networks_info,trunc_network_info,y_data)
    model.flags.training = false
    rng = MersenneTwister(1234)

    ps,st = Lux.setup(rng,model)

    
    
    input = [trunc_in,branch_in...]

    OpenHDF5(params_file,ps)
    model(input,ps,st)

    DNetPdeProblem(nx,ny,dx,dy,x,y,X,Y,cosX,sinX,bearing,model,ps,st)
end


struct PdeProblem <: AbstractPdeProblem
    nx :: Int
    ny :: Int
    dx :: Float32
    dy :: Float32
    
    x :: Array{Float32}
    y :: Array{Float32}
    X :: Array{Float32}
    Y :: Array{Float32}

    cosX :: Array{Float32}
    sinX :: Array{Float32}
    
    bearing  :: HD_Bearing
end

function PdeProblem(nx::Int,ny::Int,bearing::HD_Bearing)
    dx :: Float32 = 2pi/(nx-1)
    #dy :: Float32 = bearing.B/(ny-1)
    dy :: Float32 = 1/(ny-1)

    x = collect(LinRange(0,2pi,nx)) .|> Float32
    y = collect(LinRange(-1,1,ny)) .|> Float32
    
    X = [xv for xv in x, _ in y]
    Y = [yv for _ in x, yv in y]

    cosX = cos.(X); sinX = sin.(X)

    PdeProblem(nx,ny,dx,dy,x,y,X,Y,cosX,sinX,bearing)
end

