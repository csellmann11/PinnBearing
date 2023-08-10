function integrate_pressure(prob::DNetPdeProblem,
    pressure::Array{T}, p_fak, Φ::T) where T <: Real
    
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

    return fx,fy,Mx,My
end