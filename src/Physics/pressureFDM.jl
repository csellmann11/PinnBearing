using SparseArrays, Trapz
include("fill_mat.jl")

function spaltfunc(state_vec,pde_prob)
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
    dHdY::Matrix{Float64} = @. eps_ * pde_prob.Y .* cos(Θ - alpha_WL) # 1/2 is missing because the true y is in (0,1)
    d2HdX2::Matrix{Float64} = @. - E * cos.(Θ) - 1/2 * eps_ * pde_prob.Y .* cos(Θ - alpha_WL)

    HD ::Matrix{Float64} = @. E_ * cos.(Θ) + E * phi_ * sin.(Θ) + 1/2 * eps_dot * pde_prob.Y .* cos(Θ - alpha_WL) + 1/2 * eps_ * pde_prob.Y .* alpha_WL_dot .* sin.(Θ - alpha_WL)

    HD = HD * bearing.d/2 * 1/um
    return H[1:end-1,:],HD[1:end-1,:],dHdX[1:end-1,:],dHdY[1:end-1,:],d2HdX2[1:end-1,:]

end
function bearing_pressure(state_vec,pde_prob,dec = nothing)

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
    println("bearing width = $(bearing.B)")

    um = D/4 * om

    H,HD,dHdX,dHdY,d2Hdx2 = spaltfunc(state_vec,pde_prob)
    p_fak = C^2/(eta * abs(um) * D/2)

    anz_u = nx * (ny-2);
    nz    = 8 * nx + (ny-4) * nx * 5;

    val = Vector{T}(undef,nz) ; row = Vector{Int64}(undef,nz); col = similar(row);

    # rechte Seite 
    rhs = zeros(T,anz_u);
    
    fillMatrix!(val,row,col, ny, nx, rhs, dx, dy, H, dHdX,dHdY,d2Hdx2, HD, sign(um));
    A = sparse(row,col,val);

    if dec === nothing
        dec = lu(A)
    else
        dec = lu!(dec,A)
    end
    p_vec = dec\rhs
  
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
    fx = trapz((pde_prob.x,pde_prob.y),cos.(X) .* P)/p_fak * (D/2)^2 *bearing.B/2;
    fy = trapz((pde_prob.x,pde_prob.y),sin.(X) .* P)/p_fak * (D/2)^2 *bearing.B/2;

    return [fx,fy],dec

end