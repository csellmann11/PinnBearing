
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

function process_input(prob::DNetPdeProblem,state_vector)
    bearing = prob.bearing

    ω = bearing.om 
    
    xs,ys,dxs,dys,alpha_WL,eps_,eps_dot,alpha_WL_dot = state_vector 
    u_m = bearing.rI * ω/2

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
    
    in1,in2,in3,in4,in5,in6 = [2*E-1], [2*sqrt(_E^2 + (E * (_Φ-1))^2)/λ - 1], ϕ01_in,[2*D_m - 1],alpha_WL_in,ϕ23_in

    λ == zero(typeof(λ)) ? in2 = [0.0f0] : nothing
    width = bearing.B; width_in = [2/3 * (width - 2.5)] .|> Float32

    inputs::Vector = [in1,in1,in2,in3,in4,in5,in6, width_in]


    p_fak = λ* bearing.eta * u_m * bearing.rI^2/(bearing.c^2) * bearing.b/2 

    return inputs, E, eps_ , Φ, p_fak, alpha_WL
end


"""
    forces_dl(prob::PDE_Prob,state_vector)

Calculate the forces acting on the bearing.

# Arguments
- `prob::PDE_Prob`: The problem definition
- `state_vector`: The state vector
    contains: xs[m],ys[m],dxs[m/s],dys[m/s],alpha_WL[rad],eps_[-],eps_dot[-/s],alpha_WL_dot[-/s]

# Returns
- `fx`: The force in x direction
- `fy`: The force in y direction
"""
function forces_dl(prob::DNetPdeProblem,state_vector; Benchmark = false, pressure_return = false, parallel = false)
    model, ps, st = prob.model, prob.model_pars, prob. model_state

    if parallel 
        null = zero(eltype(state_vector))
        state_vector = [state_vector...,null,null,null+eps(Float32),null+eps(Float32)]

    end

    if Benchmark == false
         inputs, E,eps_, Φ, p_fak , alpha_WL = process_input(prob,state_vector)
    else
        inputs, E,eps_, Φ, p_fak , alpha_WL = _process_input(prob,state_vector)
    end

    net_out = model(inputs,ps,st)

    H = @. 1 + E*prob.cosX 
    
    if eps_ != zero(typeof(eps_))
        H += eps_ * 1/2 * prob.Y * cos(prob.X - alpha_WL)
    end
    pressure = reshape(net_out,prob.nx,prob.ny)./H.^2	

    pressure[pressure .< 0] .= 0

    lever = prob.Y * prob.bearing.b/2

    p_cos_x = pressure .* prob.cosX
    p_sin_x = pressure .* prob.sinX
    fx = trapz((prob.x,prob.y),p_cos_x) * p_fak
    fy = trapz((prob.x,prob.y),p_sin_x) * p_fak
    My = trapz((prob.x,prob.y),p_cos_x .* lever) * p_fak
    Mx = trapz((prob.x,prob.y),p_sin_x .* lever) * p_fak

    Rot = [cos(Φ) -sin(Φ); sin(Φ) cos(Φ)]
    F = Rot * [fx;fy]
    M = Rot * [Mx;My]
    fx,fy = F[1],F[2]
    Mx,My = M[1],M[2]

    if pressure_return
        return fx,fy,Mx,My,prob.pressure
    end
    return [fx,fy,Mx,My]
    
end