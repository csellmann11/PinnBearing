function process_input(prob::PDE_Prob,state_vector)
    bearing = prob.bearing

    ω = bearing.om |> Float32
    
    xs,ys,dxs,dys,alpha_WL,eps_,eps_dot,alpha_WL_dot = state_vector .|> Float32
    u_m = bearing.rI * ω/2

    E ::Float32 = sqrt(xs^2 + ys^2)/bearing.c
    Φ ::Float32 = atan(ys,xs) 

    _ϕ ::Float32 = (dys*xs - dxs*ys) / (xs^2 + ys^2)
    _e ::Float32 = (dxs*xs + dys*ys) / sqrt(xs^2 + ys^2)

    _Φ ::Float32 = _ϕ * bearing.rI/u_m
    _E ::Float32 = _e * bearing.rI/u_m/bearing.c 

    Eps_dot = eps_dot * bearing.rI/u_m
    Alpha_WL_dot = alpha_WL_dot * bearing.rI/u_m

    eps_max :: Float32 = 2 * (sqrt(1 - E^2*sin(alpha_WL)^2) - E * abs(cos(alpha_WL)))
    D_m = eps_/eps_max

    ϕ01 = atan((E * (_Φ-1)),_E)

    ## Schiefstellung

    f2 = Alpha_WL_dot * eps_ - eps_
    f3 = Eps_dot  
    ϕ23 = atan(f2,f3)
    #######################

    alpha_WL_in = [cos(alpha_WL),sin(alpha_WL)]
    ϕ01_in :: Vector{Float32} = [cos(ϕ01),sin(ϕ01)]
    ϕ23_in :: Vector{Float32} = [cos(ϕ23),sin(ϕ23)]

    λ::Float32 = sqrt(_E^2 + (E * (_Φ-1))^2) + sqrt(f2^2 + f3^2)
    
    in1,in2,in3,in4,in5,in6 = [2*E-1], [2*sqrt(_E^2 + (E * (_Φ-1))^2)/λ - 1], ϕ01_in,[2*D_m - 1],alpha_WL_in,ϕ23_in
    width = bearing.B; width_in = [2/3 * (width - 2.5)] .|> Float32

    inputs::Vector{Union{Vector{Float32},Float32}}   = [in1,in1,in2,in3,in4,in5,in6, width_in]

    p_fak = λ* bearing.eta * u_m * bearing.rI^2/
        (bearing.c^2) * bearing.b/2

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
function forces_dl(prob::PDE_Prob,state_vector)
    model, ps, st = prob.model, prob.model_pars, prob. model_state

    inputs, E,eps_, Φ, p_fak , alpha_WL = process_input(prob,state_vector)

    model(inputs,ps,st)

    @. prob.H = 1 + E*prob.cosX + eps_ * 1/2 * prob.Y * cos(prob.X - alpha_WL)
    
    prob.pressure .= reshape(model.output,prob.nx,prob.ny)

        
    @. prob.pressure = max(prob.pressure/prob.H^2,0) 

    fx = trapz((prob.x,prob.y),prob.pressure .* prob.cosX) * p_fak
    fy = trapz((prob.x,prob.y),prob.pressure .* prob.sinX) * p_fak

    Rot = [cos(Φ) -sin(Φ); sin(Φ) cos(Φ)]
    F = Rot * [fx;fy]
    fx,fy = F[1],F[2]

    return fx,fy
end