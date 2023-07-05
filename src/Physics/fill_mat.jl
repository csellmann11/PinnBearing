function fillMatrix!(val,row,col, ny, nx, F, dx, dy, H, dHdX, dHdY, d2Hdx2, HD, um)


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
          val[index] = h_c/dy^2 - dhdy_c/(2*dy);
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
        val[index] = -2 * (h_c * (1/dx^2 + 1/dy^2) + d2hdx2_c);
        row[index] = k;
        col[index] = l;
      
        F[k] += 12 * um * dhdx_c;

        F[k] += 12 * HD[i,j];

        ################################################
        # Sternpunkt west P_(i,j-1)
        ################################################
        is = i;

        js = j-1;
       
        l = (is) + (nx) * (js - 2);

  
        if (js > 1) 
          index += 1;
          val[index] = h_c/dy^2 + dhdy_c/(2*dy);
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