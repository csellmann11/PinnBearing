"""
    fillMatrix!(val, row, col, ny, nx, F, dx, dy, H, dHdX, dHdY, d2Hdx2, rhs_pde, alpha)

Fill the matrix values, rows, and columns for a finite difference discretization of a PDE.

# Arguments
- `val`: Array to store the non-zero values of the matrix.
- `row`: Array to store the row indices of the non-zero values.
- `col`: Array to store the column indices of the non-zero values.
- `ny`: Number of grid points in the y-direction.
- `nx`: Number of grid points in the x-direction.
- `F`: Right-hand side of the PDE.
- `dx`: Grid spacing in the x-direction.
- `dy`: Grid spacing in the y-direction.
- `H`: Array representing the clearance hight.
- `dHdX`: Array representing the x-derivative of the field `H`.
- `dHdY`: Array representing the y-derivative of the field `H`.
- `d2Hdx2`: Array representing the second x-derivative of the field `H`.
- `rhs_pde`: Right-hand side values of the PDE at each grid point.
- `alpha`: (rI/B)^2, scaled bearing radius squared divided by the bearing width squared.

# Description
This function fills the `val`, `row`, and `col` arrays based on the finite difference discretization 
of a given PDE. The function uses a loop over the grid points and computes the matrix entries 
based on the provided field values and their derivatives.

# Returns
The function modifies the `val`, `row`, `col`, and `F` arrays in-place.
"""
function fillMatrix!(val,row,col, ny, nx, F, dx, dy, H, dHdX, dHdY, d2Hdx2, rhs_pde, alpha)
    index = 0;
    for j ∈ 2:(ny-1), i ∈ 1:nx

        im1 = i - 1;
        if (im1 < 1) 
          im1 = nx;
        end

        ip1 = i + 1;
        if (ip1 > nx) 
          ip1 = 1;
        end

        h_c         = H[i,j];
        dhdx_c      = dHdX[i,j];
        dhdy_c      = dHdY[i,j];
        d2hdx2_c    = d2Hdx2[i,j];

        k = (i) + (nx) * (j-2);

        ################################################
        # Sternpunkt north P_(i-1,j)
        ################################################
        is = i-1;
        if (is < 1) 
            is = nx;
        end
        js = j;

        l = (is) + (nx) * (js - 2);

        index += 1;

        val[index] = h_c/dx^2 + dhdx_c/(2*dx); 
        row[index] = k; 
        col[index] = l;

        ################################################
        # Sternpunkt east P_(i,j+1)
        ################################################
        is = i;

        js = j + 1;

        
        l = (is) + (nx) * (js - 2);

        
        if ((j + 1) < ny) 
          
          index += 1;
          val[index] = alpha * (h_c/dy^2 - dhdy_c/(2*dy));
          row[index] = k; 
          col[index] = l;
        end
        ################################################
        # Sternpunkt mitte
        ################################################
        is = i;

        js = j;


        l = (is) + (nx) * (js - 2);

        index += 1;
        val[index] = -2 * (h_c * (1/dx^2 + alpha/dy^2) + d2hdx2_c);
        row[index] = k;
        col[index] = l;
      
        F[k] = rhs_pde[i,j];
        ################################################
        # Sternpunkt west P_(i,j-1)
        ################################################
        is = i;

        js = j-1;
       
        l = (is) + (nx) * (js - 2);

  
        if (js > 1) 
          index += 1;
          val[index] = alpha * (h_c/dy^2 + dhdy_c/(2*dy));
          row[index] = k;
          col[index] = l;
        end

        ################################################
        # Sternpunkt south
        ################################################
        is = i+1;
 
        if (is > nx) 
          is = 1;
        end

        js = j;

        l = (is) + (nx) * (js - 2);

        index += 1;
        val[index] = h_c/dx^2 - dhdx_c/(2*dx);
        row[index] = k;
        col[index] = l;

    end
end