
"""
    _process_input(prob::PDE_Prob,state_vector)

# Warning: 
- This function is only for testing purposes.
Use `process_input` instead.

"""

function process_input_nl(prob::DNetPdeProblem,state_vector, misa)
    bearing = prob.bearing

    ω = bearing.om |> Float32
    
    if misa
        xs,ys,alpha_WL,eps_ = state_vector .|> Float32
    else
        xs,ys = state_vector .|> Float32
    end

    omega_max = bearing.omega_max |> Float32
    E ::Float32 = sqrt(xs^2 + ys^2)/bearing.c
    #println("E: ",E)
    Φ ::Float32 = atan(ys,xs) 

    if misa
        alpha_WL_in = [cos(alpha_WL),sin(alpha_WL)]

        in1,in2,in3,in4 = [2*E-1], [2ω/omega_max-1], [eps_], [alpha_WL_in]
        inputs  = [in1,in1,in2,in3,in4]
    else
        alpha_WL = 0.0
        eps_ = 0.0
        in1,in2 = [2*E-1.0f0], [2*ω/omega_max-1.0f0]
        inputs   = [in1,in1,in2]
    end

    p_fak = bearing.b/2 * bearing.rI * bearing.p_atm

    return inputs, E, eps_ , Φ, alpha_WL, p_fak
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
function forces_dl_nl(prob::DNetPdeProblem,state_vector; pressure_return = false, misa = false)
    model, ps, st = prob.model, prob.model_pars, prob. model_state

   
    inputs, E,eps_, Φ, alpha_WL, p_fak = process_input_nl(prob,state_vector, misa)
    

    model(inputs,ps,st)

    #@. prob.H = 1 + E*prob.cosX + eps_ * 1/2 * prob.Y * cos(prob.X - alpha_WL)
    
    prob.pressure .= reshape(model.output,prob.nx,prob.ny) .+ 1

    fx = trapz((prob.x,prob.y),prob.pressure .* prob.cosX) * p_fak
    fy = trapz((prob.x,prob.y),prob.pressure .* prob.sinX) * p_fak

    Rot = [cos(Φ) -sin(Φ); sin(Φ) cos(Φ)]
    F = Rot * [fx;fy]
    fx,fy = F[1],F[2]

    if pressure_return
        return fx,fy,prob.pressure
    end
    return fx,fy
end