"""
    _process_input(prob::PDE_Prob,state_vector)

# Warning: 
- This function is only for testing purposes.
Use `process_input` instead.

"""
function _process_input(prob::DNetPdeProblem,state_vector)
    bearing = prob.bearing

    ω = bearing.om 
    
    E,A1,ϕ01,D_m,alpha_WL,ϕ23 = state_vector .|> Float32
    u_m = bearing.rI * ω/2
    Φ = 0.0

    eps_max = 2 * (sqrt(1 - E^2 * sin(alpha_WL)^2) - E * abs(cos(alpha_WL)))
    eps_ = D_m * eps_max

    #######################

    alpha_WL_in :: Vector{Float32}  = [cos(alpha_WL),sin(alpha_WL)]
    ϕ01_in :: Vector{Float32}       = [cos(ϕ01),sin(ϕ01)]
    ϕ23_in :: Vector{Float32}       = [cos(ϕ23),sin(ϕ23)]

    λ::Float32 = 1.0
    
    in1,in2,in3,in4,in5,in6 = [2*E-1], [2*A1 - 1], ϕ01_in,[2*D_m - 1],alpha_WL_in,ϕ23_in
    width = bearing.B; width_in = [2/3 * (width - 2.5)] .|> Float32

    inputs::Vector{Union{Vector{Float32},Float32,Nothing}} = [nothing,in1,in2,in3,in4,in5,in6, width_in]

    p_fak = λ* bearing.eta * u_m * bearing.rI^2/
        (bearing.c^2) * bearing.b/2

    return inputs, E, eps_ , Φ, p_fak, alpha_WL
end

"""
process_input(prob::DNetPdeProblem, 
    state_vector::Vector{T}, benchmark::Bool) where T <: Real

# Arguments
- `prob::PDE_Prob`: The problem definition
- `state_vector`: The state vector
    contains: xs[m],ys[m],dxs[m/s],dys[m/s],alpha_WL[rad],eps_[-],eps_dot[-/s],alpha_WL_dot[-/s]
- `Benchmark`: If true, the function `_process_input` is used instead of `process_input` (keep it false for use)
# Returns
- `fx`: The force in x direction
- `fy`: The force in y direction
# Description
This function processes the state vector and returns the inputs for the neural network.

# Returns
- `inputs`: The inputs for the neural network
- `E`: The eccentricity ratio
- `eps_`: First factor influencing the rotor's tilt
- `Φ`: The angle between the x-axis and the eccentricity vector
- `p_fak`: The pressure factor
- `alpha_WL`: second factor influencing the rotor's tilt
"""
function process_input(prob::DNetPdeProblem, 
    state_vector::Vector{T}, benchmark::Bool) where T <: Real

    if benchmark
        return _process_input_bench(prob,state_vector)
    end

    length(state_vector) == 4 ? parallel = true : false

    if parallel # if parallel add zeros to the state vector
        null = zero(eltype(state_vector))
        state_vector = [state_vector...,null,null,null+eps(Float32),null+eps(Float32)]
    end
    @assert length(state_vector) == 8 "The length of the state vector is not 8 or 4"

    bearing = prob.bearing
    
    xs,ys,dxs,dys,alpha_WL,eps_,eps_dot,alpha_WL_dot = state_vector 
    u_m = bearing.rI * bearing.om/2

    # Make inputs dimensionless
    E = sqrt(xs^2 + ys^2)/bearing.c
    Φ = atan(ys,xs) 
    _ϕ = (dys*xs - dxs*ys) / (xs^2 + ys^2)
    _e = (dxs*xs + dys*ys) / sqrt(xs^2 + ys^2)
    _Φ = _ϕ * bearing.rI/u_m
    _E = _e * bearing.rI/u_m/bearing.c 

    Eps_dot = eps_dot * bearing.rI/u_m
    Alpha_WL_dot = alpha_WL_dot * bearing.rI/u_m

    eps_max = 2 * (sqrt(1 - E^2*sin(alpha_WL)^2) - E * abs(cos(alpha_WL)))
    eps_max != zero(typeof(eps_max)) ? D_m = eps_/eps_max : D_m = 0.0f0 
    ϕ01 = atan((E * (_Φ-1)),_E)

    ## Schiefstellung

    f2 = Alpha_WL_dot * eps_ - eps_
    f3 = Eps_dot  
    ϕ23 = atan(f2,f3)
    #######################

    alpha_WL_in = [cos(alpha_WL),sin(alpha_WL)]
    ϕ01_in  = [cos(ϕ01),sin(ϕ01)]
    ϕ23_in  = [cos(ϕ23),sin(ϕ23)]

    λ = sqrt(_E^2 + (E * (_Φ-1))^2) + sqrt(f2^2 + f3^2)
    
    in1,in2,in3,in4,in5,in6 = [2*E-1], [2*sqrt(_E^2 + (E * (_Φ-1))^2)/λ - 1], ϕ01_in,[2*D_m - 1],
                    alpha_WL_in,ϕ23_in

    λ == zero(typeof(λ)) ? in2 = [0.0f0] : nothing
    width = bearing.B; width_in = [2/3 * (width - 2.5)] .|> Float32

    inputs::Vector = [in1,in1,in2,in3,in4,in5,in6, width_in]


    p_fak = λ* bearing.eta * u_m * bearing.rI^2/(bearing.c^2) * bearing.b/2 

    return inputs, E, eps_ , Φ, p_fak, alpha_WL
end

function pos(x::T) where T <: Real
    null = zero(T)
    if x > null
        return x
    else
        return null
    end
end

"""
    forces_dl(prob::PDE_Prob,state_vector)

Calculate the forces acting on the bearing.

# Arguments
- `prob::PDE_Prob`: The problem definition
- `state_vector`: The state vector
    contains: xs[m],ys[m],dxs[m/s],dys[m/s],alpha_WL[rad],eps_[-],eps_dot[-/s],alpha_WL_dot[-/s]
- `Benchmark`: If true, the function `_process_input` is used instead of `process_input` (keep it false for use)
- `pressure_return`: If true, the pressure is returned as well
# Returns
- `fx`: The force in x direction
- `fy`: The force in y direction
"""
function forces_dl(state_vector::Vector{T},prob::DNetPdeProblem; 
    benchmark = false, 
    pressure_return = false) where T <: Real

    model, ps, st = prob.model, prob.model_pars, prob. model_state

    inputs, E,eps_, Φ, p_fak , alpha_WL = process_input(prob,state_vector,benchmark)

    net_out = model(inputs,ps,st) # forward pass

    if length(state_vector) == 4
        H = @. 1 + E*prob.cosX 
    else
        H = @. 1 + E*prob.cosX + eps_ * 1/2 * prob.Y * cos(prob.X - alpha_WL)
    end

    # norm pressure and gumbel cavitation model
    pressure = reshape(net_out,prob.nx,prob.ny)./H.^2 .|> pos

    # calculation of the forces
    fx,fy,Mx,My = integrate_pressure(prob,pressure,p_fak,Φ)

    if pressure_return
        return fx,fy,Mx,My,pressure
    end
    return [fx,fy,Mx,My]
    
end


