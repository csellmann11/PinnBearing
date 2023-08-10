include("utils_force.jl")
"""
    _process_input(prob::PDE_Prob,state_vector)

# Warning: 
- This function is only for testing purposes.
Use `process_input` instead.

"""

function process_input_nl(prob::DNetPdeProblem,state_vector)
    bearing = prob.bearing

    ω = bearing.om |> Float32
    
    xs,ys,alpha_WL,eps_ = state_vector .|> Float32
    
    omega_max = bearing.omega_max |> Float32
    E ::Float32 = sqrt(xs^2 + ys^2)/bearing.c

    Φ ::Float32 = atan(ys,xs) 
 
    alpha_WL_in = [cos(alpha_WL),sin(alpha_WL)]

    in1,in2,in3,in4 = [2*E-1], [2ω/omega_max-1], [eps_], [alpha_WL_in]
    inputs  = [in1,in1,in2,in3,in4]
   

    p_fak = bearing.b/2 * bearing.rI * bearing.p_atm

    return inputs, Φ, p_fak
end


"""
    forces_dl(prob::PDE_Prob,state_vector)

Calculate the static forces acting on the bearing.

# Arguments
- `prob::PDE_Prob`: The problem definition
- `state_vector`: The state vector
    contains: xs[m],ys[m],alpha_WL[rad],eps_[-] or xs[m],ys[m]

# Returns
- `fx`: The force in x direction
- `fy`: The force in y direction
"""
function forces_dl_nl(state_vector::Vector{T},prob::DNetPdeProblem; pressure_return = false) where T <: Real

    model, ps, st = prob.model, prob.model_pars, prob. model_state

    l_state = length(state_vector)
    if l_state == 2
        state_vector = [state_vector...,zero(T),zero(T)]
    elseif l_state == 4
        @warn("net only trained for 2 state variables so far")
    else
        error("state_vector must have length 2 or 4")
    end

    inputs, Φ, p_fak = process_input_nl(prob,state_vector)
    

    net_out = model(inputs,ps,st)
    
    pressure = reshape(net_out,prob.nx,prob.ny) .+ 1

    fx,fy,Mx,My = integrate_pressure(prob,pressure,p_fak,Φ)

    if pressure_return
        return fx,fy,Mx,My,pressure
    end
    return fx,fy,Mx,My
end

