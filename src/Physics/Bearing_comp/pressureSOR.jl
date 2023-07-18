using LoopVectorization
using Trapz
using Statistics
function calc_press(P, P_new, H, θ, nx, ny, Δx, Δy, Ks, p0, c, e, iter, ω_relax, alpha, diff, Λ_c)
    
    e_g = e
    for i in 1:iter

        if i < 6500
           
            e = 0.5 + (e_g-0.5) * i/6500
        end
        
        for i in 1:nx-1
            for j in 2:ny-1

                A = 1/8 * (H[mod1(i+1,nx-1),j] + H[i,j])^3 * (P[mod1(i+1,nx-1),j] + P[i,j])
                B = 1/8 * (H[i,j] + H[mod1(i-1,nx-1),j])^3 * (P[i,j] + P[mod1(i-1,nx-1),j])
                C = alpha/8 * (Δx/Δy)^2 * (H[i,j+1] + H[i,j])^3 * (P[i,j+1] + P[i,j])
                D = alpha/8 * (Δx/Δy)^2 * (H[i,j] + H[i,j-1])^3 * (P[i,j] + P[i,j-1])
                E = Λ_c * Δx * (H[mod1(i+1,nx-1),j] - H[mod1(i-1,nx-1),j])
                
                s1 = -Λ_c * H[i,j] * (P[mod1(i+1,nx-1),j] - P_new[mod1(i-1,nx-1),j])* Δx
                s2 = A * P[mod1(i+1,nx-1),j] + B * P_new[mod1(i-1,nx-1),j] + C * P[i,j+1] + D * P[i,j-1]

                @inbounds P_new[i,j] = (1- ω_relax) * P[i,j] + ω_relax * (s1 + s2)/(A + B + C + D + E)



            end
            P_row = @view P_new[i,:]
            mean_pressure = mean(P_row) - 1
            @. H[i,:] = 1 .+ e .* cos.(@views θ[i,:]) .+ mean_pressure .* p0 ./ Ks #./ c

        end
        
        @. diff .= P - P_new
        diff .= abs.(diff)
   
        abs_error = maximum(diff)
        diff .= diff ./ P
        rel_error = maximum(diff)

        if abs_error < 1e-7 && rel_error < 1e-7
            break
        end

        if i == iter
            println("Did not converge")
            println("abs_error = $abs_error, rel_error = $rel_error")
            break
        end

        P .= P_new
    end
end

function nonlinear_pressure(state_vec,prob)
    bearing = prob.bearing

    η = bearing.eta; ω = bearing.om; r_I = bearing.rI; b = bearing.b
    c = bearing.c; um = 1/2 * ω * r_I; p0 = bearing.p_atm

    Λ_c = 6 * η * ω * r_I^2/(p0 * c^2); 
    T  = 6 * η * r_I^2/(p0 * c^2);
    alpha = (r_I/b)^2

    nx = 64; ny = 24; 

    xs,ys = state_vec
    e = sqrt(xs^2 + ys^2)/c
    ϕ = atan(ys,xs)

    x = LinRange(0, 2π, nx); y = LinRange(0, 1, ny)
    Δx = x[2] - x[1]; Δy = y[2] - y[1]

    X = zeros(nx-1, ny); Y = zeros(nx-1, ny)
    for (i,xv) in enumerate(x[1:end-1])
        for (j,yv) in enumerate(y)
            X[i, j] = xv
            Y[i, j] = yv
        end
    end
    

    θ = X .- ϕ 
    H = 1 .+ 0.5 * cos.(θ)  
    P = ones(nx-1,ny)
    P_new = ones(axes(P))
    diff = zeros(axes(P))


    Ks = bearing.Ksc
    

    ω_relax = 0.1
    iter = 30000

    calc_press(P, P_new, H, θ, nx, ny, Δx, Δy, Ks, p0, c, e, iter, ω_relax, alpha, diff, Λ_c)

    pressure = [P; reshape(P[end,:],1,:)]

    X_full = [X; reshape(X[end,:],1,:)]

    fx = trapz((x,y),pressure .* cos.(X_full))*b * r_I * p0; 
    fy = trapz((x,y),pressure .* sin.(X_full))*b * r_I * p0;

    return fx, fy
end

