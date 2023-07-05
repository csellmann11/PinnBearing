using Lux

########################
# Structure Definition #
########################
struct FFCN{L<:NamedTuple,T<:NamedTuple} <: Lux.AbstractExplicitContainerLayer{(:layers,:pre_layers)}
    layers::L
    pre_layers::T
end

function FFCN(in_dim, out_dim, n_hidden_layers, size_hidden, activation = tanh)
    
    n_layers = n_hidden_layers + 2
    n_pre_layers = 2


    layers = Array{Pair}(undef, n_layers)
    layers[1] = Pair(Symbol("layer_",1),Dense(in_dim, size_hidden, activation))
    for i in 1:n_hidden_layers
        layers[i+1] = Pair(Symbol("layer_",i+1),Dense(size_hidden, size_hidden, activation))
    end
    layers[n_layers] = Pair(Symbol("layer_",n_layers),Dense(size_hidden, out_dim, identity))

    pre_layers = Array{Pair}(undef, n_pre_layers)
    pre_layers[1] = Pair(Symbol("layer_",1),Dense(in_dim, size_hidden, activation))
    pre_layers[2] = Pair(Symbol("layer_",2),Dense(in_dim, size_hidden, activation))

    return FFCN(NamedTuple(layers), NamedTuple(pre_layers))
end

######################
# Forward Evaluation #
######################

function (net::FFCN)(x, ps::NamedTuple,st)


    u = net.pre_layers.layer_1(x, ps.pre_layers.layer_1, st.pre_layers.layer_1)[1]
    v = net.pre_layers.layer_2(x, ps.pre_layers.layer_2, st.pre_layers.layer_2)[1]
    

    pre_faktor = 1 .- u .+ v

    for (idx,layer) in enumerate(net.layers)
        if idx != length(net.layers)
            x = layer(x, ps.layers[idx], st.layers[idx])[1]
            x = x .* pre_faktor
            
        else
            x = layer(x, ps.layers[idx], st.layers[idx])[1]
        end
    end

    return x, st
end 




