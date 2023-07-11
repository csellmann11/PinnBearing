using SparseArrays, Trapz
include("fill_mat.jl")

"""
spaltfunc(state_vec,pde_prob::AbstractPdeProblem)

Calculate the clearance function and its derivatives.

# Arguments
- `state_vector`: The state vector
    contains: xs[m],ys[m],dxs[m/s],dys[m/s],alpha_WL[rad],eps_[-],eps_dot[-/s],alpha_WL_dot[-/s]
- `pde_prob::AbstractPdeProblem`: The problem definition

# Returns
- `H`: The clearance function
- `HD`: The derivative of the clearance function with respect to time
- `dHdX`: The derivative of the clearance function with respect to the x-coordinate
- `dHdY`: The derivative of the clearance function with respect to the y-coordinate
- `d2HdX2`: The second derivative of the clearance function with respect to the x-coordinate
"""
function spaltfunc(state_vec,pde_prob::AbstractPdeProblem)
    bearing = pde_prob.bearing
    um = bearing.d/4 * bearing.om
    x,y,dx,dy,alpha_WL,eps_,eps_dot,alpha_WL_dot = state_vec .|> Float64

    E = sqrt(x^2 + y^2)/bearing.c
    phi = atan(y,x)
    phi_ = (dy*x -y * dx) / (x^2 + y^2)
    E_   = (dy*y + x * dx) / sqrt(x^2 + y^2)/ bearing.c

    Θ = pde_prob.X .- phi
    H  ::Matrix{Float64}  = @. 1 .+ E * cos.(Θ) + 1/2 * eps_ * pde_prob.Y .* cos(Θ - alpha_WL)
    dHdX::Matrix{Float64} = @. - E * sin.(Θ)  - 1/2 * eps_ * pde_prob.Y .* sin.(Θ - alpha_WL)
    dHdY::Matrix{Float64} = @. eps_ * cos(Θ - alpha_WL)# 1/2 is missing because the true y is in (0,B)
    d2HdX2::Matrix{Float64} = @. - E * cos.(Θ) - 1/2 * eps_ * pde_prob.Y .* cos(Θ - alpha_WL)

    HD ::Matrix{Float64} = @. E_ * cos.(Θ) + E * phi_ * sin.(Θ) + 1/2 * eps_dot * pde_prob.Y .* cos(Θ - alpha_WL) + 1/2 * eps_ * pde_prob.Y .* alpha_WL_dot .* sin.(Θ - alpha_WL)

    HD = HD * bearing.d/2 * 1/um
    return H[1:end-1,:],HD[1:end-1,:],dHdX[1:end-1,:],dHdY[1:end-1,:],d2HdX2[1:end-1,:]

end

"""
bearing_pressure(state_vec,pde_prob::AbstractPdeProblem)

Calculate the forces in the bearing with FDM.

# Arguments
- `state_vector`: The state vector
    contains: xs[m],ys[m],dxs[m/s],dys[m/s],alpha_WL[rad],eps_[-],eps_dot[-/s],alpha_WL_dot[-/s]
- `pde_prob::AbstractPdeProblem`: The problem definition

# Returns
- `fx`: The force in x-direction
- `fy`: The force in y-direction
"""
function bearing_pressure(state_vec,pde_prob::AbstractPdeProblem; pressure_ret = false)

    T = Float64
    bearing = pde_prob.bearing

    D   =bearing.d;
    B   =bearing.b;
    C   =bearing.c;
    eta =bearing.eta;
    om  =bearing.om;

    nx = pde_prob.nx-1
    ny = pde_prob.ny
    dx ::Float64 = pde_prob.dx
    dy ::Float64 = pde_prob.dy

    um = D/4 * om

    H,HD,dHdX,dHdY,d2Hdx2 = spaltfunc(state_vec,pde_prob)
    p_fak = C^2/(eta * abs(um) * D/2)

    anz_u = nx * (ny-2);
    nz    = 8 * nx + (ny-4) * nx * 5;

    val = Vector{T}(undef,nz) ; row = Vector{Int64}(undef,nz); col = similar(row);

    # rechte Seite 
    rhs = zeros(T,anz_u);
    
    rhs_pde = 12 * (sign(um) * dHdX + HD)
    fillMatrix!(val,row,col, ny, nx, rhs, dx, dy, H, dHdX,dHdY,d2Hdx2, rhs_pde, bearing.alpha);
    A = sparse(row,col,val);

    # p_vec = klu(A)\rhs
    # p_vec = p_vec ./ reshape(H[:,2:end-1],:,1).^2

    if bearing.sysMat_dec === nothing
        bearing.sysMat_dec = klu(A)
        println("Decomposition done")
    end

    p_vec = klu!(bearing.sysMat_dec,A)\rhs
    p_vec = p_vec ./ reshape(H[:,2:end-1],:,1).^2


    P = zeros(T,nx+1,ny)
    for i ∈ eachindex(p_vec)
        if p_vec[i] < 0
            p_vec[i] = 0
        end
    end
    for i ∈ 1:nx, j ∈ 1:ny
        if j == 1
            P[i,j] = 0;
        elseif j == ny
            P[i,j] = 0;
        else
            P[i,j] = p_vec[(j-2)*nx + i];
        end
    end

    P[nx+1,:] = P[1,:];

    X = pde_prob.X
    fx = trapz((pde_prob.x,pde_prob.y),cos.(X) .* P)/p_fak * (D/2) * bearing.b/2;
    fy = trapz((pde_prob.x,pde_prob.y),sin.(X) .* P)/p_fak * (D/2) * bearing.b/2;

    if pressure_ret
        return [fx,fy],P
    end

    return [fx,fy]

end




#######################################
# Version with DL inputs
#######################################

"""
_spaltfunc(state_vec,pde_prob::AbstractPdeProblem)

# Warning: 
- This function is not intended to be used directly. Use `spaltfunc` instead.
    Only for testing purposes.
"""
function _spaltfunc(state_vec,pde_prob::AbstractPdeProblem)
    E,A1,ϕ01,Dm,alpha_WL,ϕ23 = state_vec .|> Float64

    A2 = 1 - A1

    eps_max = 2 * (sqrt(1 - E^2 * sin(alpha_WL)^2) - E * abs(cos(alpha_WL)))
    eps_ = Dm * eps_max
    phi = 0.0

    Θ = pde_prob.X .- phi

    H  ::Matrix{Float64}    = @. 1 .+ E * cos.(Θ) + 1/2 * eps_ * pde_prob.Y .* cos(Θ - alpha_WL)
    dHdX::Matrix{Float64}   = @.    - E * sin.(Θ) - 1/2 * eps_ * pde_prob.Y .* sin.(Θ - alpha_WL)
    d2HdX2::Matrix{Float64} = @.    - E * cos.(Θ) - 1/2 * eps_ * pde_prob.Y .* cos(Θ - alpha_WL)

    dHdY::Matrix{Float64} = @. eps_ * cos(Θ - alpha_WL)  # 1/2 is missing because the true y is in (0,bearing.B)

    pde_rhs = @. 12 * (A1 * cos.(Θ .- ϕ01) + 1/2 * A2 * pde_prob.Y .* cos.(Θ .- ϕ23 .- alpha_WL))
    
    return H[1:end-1,:],pde_rhs[1:end-1,:],dHdX[1:end-1,:],dHdY[1:end-1,:],d2HdX2[1:end-1,:]

end

"""
spaltfunc(state_vec,pde_prob::AbstractPdeProblem)

# Warning:
- This function is not intended to be used directly. Use `bearing_pressure` instead.
    Only for testing purposes.
"""
function _bearing_pressure(state_vec,pde_prob)

    T = Float64
    bearing = pde_prob.bearing

    D   =bearing.d;
    B   =bearing.b;
    C   =bearing.c;
    eta =bearing.eta;
    om  =bearing.om;

    nx = pde_prob.nx-1
    ny = pde_prob.ny
    dx ::Float64 = pde_prob.dx
    dy ::Float64 = pde_prob.dy

    um = D/4 * om

    H,pde_rhs,dHdX,dHdY,d2Hdx2 = _spaltfunc(state_vec,pde_prob)
    p_fak = C^2/(eta * abs(um) * D/2)

    anz_u = nx * (ny-2);
    nz    = 8 * nx + (ny-4) * nx * 5;

    val = Vector{T}(undef,nz) ; row = Vector{Int64}(undef,nz); col = similar(row);

    # rechte Seite 
    rhs = zeros(T,anz_u);
    
    fillMatrix!(val,row,col, ny, nx, rhs, dx, dy, H, dHdX,dHdY,d2Hdx2, pde_rhs, bearing.alpha);
    A = sparse(row,col,val);

    

    p_vec = klu(A)\rhs
    p_vec = p_vec ./ reshape(H[:,2:end-1],:,1).^2

    P = zeros(T,nx+1,ny)
    for i ∈ eachindex(p_vec)
        if p_vec[i] < 0
            p_vec[i] = 0
        end
    end
    for i ∈ 1:nx, j ∈ 1:ny
        if j == 1
            P[i,j] = 0;
        elseif j == ny
            P[i,j] = 0;
        else
            P[i,j] = p_vec[(j-2)*nx + i];
        end
    end

    P[nx+1,:] = P[1,:];

    fx = trapz((pde_prob.x,pde_prob.y),pde_prob.cosX .* P)/p_fak * (D/2) * bearing.b/2;
    fy = trapz((pde_prob.x,pde_prob.y),pde_prob.sinX .* P)/p_fak * (D/2) * bearing.b/2;
    return [fx,fy]

end