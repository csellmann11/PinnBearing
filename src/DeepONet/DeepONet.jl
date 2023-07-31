using Lux, Random

include("FFCN.jl")

function mymul!(out::Vector{Float32},mat::Matrix{Float32},vec::Vector{Float32})
    """
    Multiplication of the transposed a matrix with a vector and store the result in a vector.
    u[i] = A[j,i] * v[j]
    """
    @assert size(mat,1) == size(vec,1) "First dimension of mat and vec must match"
    @turbo for i in axes(mat,2)
        out[i] = zero(eltype(out))	
        for j in axes(mat,1)
            out[i] += mat[j,i] * vec[j]
        end
    end
end

function mymul!(out::Vector{T},mat::Matrix{Float32},vec::Vector{T}) where {T <: Real}
    """
    Multiplication of the transposed a matrix with a vector and store the result in a vector.
    u[i] = A[j,i] * v[j]
    """
    mul!(out,mat',vec)
end

########################
# Structure Definition #
########################

mutable struct Flags
    first_call::Bool
    training::Bool
end
function Flags()
    return Flags(true,true)
end

struct DeepONet{B<:NamedTuple,T<:FFCN} <: Lux.AbstractExplicitContainerLayer{(:branch_networks,:trunc_network)}
    branch_networks::B
    trunc_network::T

    trunc_output::Array{Float32}
    distance::Array{Float32}

    num_branches::Int

    flags ::Flags
end

function DeepONet(branch_networks_info,trunc_network_info,y_data)
    
    n_branches = length(branch_networks_info)
    
    branch_nets = Array{Pair}(undef, n_branches)
    for i in eachindex(branch_nets)
        branch_nets[i] = Symbol("branch_network",i) => FFCN(branch_networks_info[i]...)
    end

    trunc_net = FFCN(trunc_network_info...)

    distance = cos.(pi/2 * y_data[:]) .|> Float32


    reduced_dim = trunc_network_info[2]; len_data = length(y_data)
    trunc_output = Array{Float32}(undef, reduced_dim, len_data)

    return DeepONet(NamedTuple(branch_nets),trunc_net,trunc_output,distance,n_branches, Flags())
end

######################
# Forward Evaluation #
######################
function DeepONet_AnalyseEval(net,x,ps::NamedTuple,st)
    trunc_output = net.trunc_network(x[1],ps.trunc_network,st.trunc_network)[1] # red_dim, n_samples

    branch_outputs = [net.branch_networks[i](x[i+1],ps.branch_networks[i],st.branch_networks[i])[1] 
                        for i in 1:net.num_branches]
    
    branch_out = reduce(.*,branch_outputs) # red_dim, batch
    return trunc_output' * branch_out, branch_outputs
end


function DeepONet_TrainEval(net,x,ps::NamedTuple,st)
    trunc_output = net.trunc_network(x[1],ps.trunc_network,st.trunc_network)[1] # red_dim, n_samples

    branch_outputs = [net.branch_networks[i](x[i+1],ps.branch_networks[i],st.branch_networks[i])[1] 
                        for i in 1:net.num_branches]
    
    branch_out = reduce(.*,branch_outputs) # red_dim, batch
    return trunc_output' * branch_out
end

function DeepONet_InferenceEval(net,x,ps::NamedTuple,st)
    if net.flags.first_call
        net.flags.first_call = false
        net.trunc_output .= net.trunc_network(x[1],ps.trunc_network,st.trunc_network)[1] .* net.distance' # red_dim, n_samples
    end

    branch_outputs = [net.branch_networks[i](x[i+1],ps.branch_networks[i],st.branch_networks[i])[1] for i in 1:net.num_branches]
    branch_out = reduce(.*,branch_outputs) # red_dim, batch

    T = eltype(branch_out)
    output = Array{T}(undef, size(net.trunc_output,2))
    mymul!(output,net.trunc_output,branch_out)

    return output
end

function (net::DeepONet)(x, ps::NamedTuple,st)
    if net.flags.training
        return DeepONet_Train(net,x,ps,st)
    else
        output = DeepONet_InferenceEval(net,x,ps,st) 
        return output
    end
end





