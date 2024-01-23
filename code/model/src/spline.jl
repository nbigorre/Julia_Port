function spline(n, x, y, b, c, d)
      #---------------------------------------------                                                                                    
      #  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed      
      #  for a cubic interpolating spline                                     
      #                                                                       
      #    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3  
      #                                                                       
      #    for  x(i) .le. x .le. x(i+1)                                       
      #                                                                       
      #  input..                                                              
      #                                                                       
      #    n = the number of data points or knots (n.ge.2)                    
      #    x = the abscissas of the knots in strictly increasing order        
      #    y = the ordinates of the knots                                     
      #                                                                       
      #  output..                                                             
      #                                                                       
      #    b, c, d  = arrays of spline coefficients as defined above.         
      #                                                                       
      #  using  p  to denote differentiation,                                 
      #                                                                       
      #    y(i) = s(x(i))                                                     
      #    b(i) = sp(x(i))                                                    
      #    c(i) = spp(x(i))/2                                                 
      #    d(i) = sppp(x(i))/6  (derivative from the right)                   
      #                                                                       
      #  the accompanying function subprogram  seval  can be used             
      #  to evaluate the spline.                                              
                                                                
      local nm1 = n - 1
      if (n < 2)
            return
      end
      if (n < 3)
            b[1] = (y[2] - y[1]) / (x[2] - x[1])
            c[1] = 0e0
            d[1] = 0e0
            b[2] = b[1]
            c[2] = 0e0
            d[2] = 0e0
            return
      end
      #                                                                       
      #  set up tridiagonal system                                            
      #                                                                       
      #  b = diagonal, d = offdiagonal, c = right hand side.                  
      #                                                                       
      d[1] = x[2] - x[1]
      c[2] = (y[2] - y[1]) / d[1]

      for i in 2:nm1
            d[i] = x[i+1] - x[i]
            b[i] = 2e0 * (d[i-1] + d[i])
            c[i+1] = (y[i+1] - y[i]) / d[i]
            c[i] = c[i+1] - c[i]
      end
      #                                                                       
      #  end conditions.  third derivatives at  x(1)  and  x(n)               
      #  obtained from divided differences                                    
      #                                                                       
      b[1] = -d[1]
      b[n] = -d[n-1]
      c[1] = 0e0
      c[n] = 0e0
      if (n != 3)
            c[1] = c[3] / (x[4] - x[2]) - c[2] / (x[3] - x[1])
            c[n] = c[n-1] / (x[n] - x[n-2]) - c[n-2] / (x[n-1] - x[n-3])
            c[1] = c[1] * d[1]^2 / (x[4] - x[1])
            c[n] = -c[n] * d[n-1]^2 / (x[n] - x[n-3])
      end
      #                                                                       
      #  forward elimination                                                  
      #                                                                       
      for i in 2:n
            local t = d[i-1] / b[i-1]
            b[i] = b[i] - t * d[i-1]
            c[i] = c[i] - t * c[i-1]
      end
      #                                                                       
      #  back substitution                                                    
      #                                                                       
      c[n] = c[n] / b[n]
      for ib in 1:nm1
            local i = n - ib
            c[i] = (c[i] - d[i] * c[i+1]) / b[i]
      end
      #                                                                       
      #  c(i) is now the sigma(i) of the text                                 
      #                                                                       
      #  compute polynomial coefficients                                      
      #                                                                       
      b[n] = (y[n] - y[nm1]) / d[nm1] + d[nm1] * (c[nm1] + 2e0 * c[n])
      for i in 1:nm1
            b[i] = (y[i+1] - y[i]) / d[i] - d[i] * (c[i+1] + 2e0 * c[i])
            d[i] = (c[i+1] - c[i]) / d[i]
            c[i] = 3e0 * c[i]
      end
      c[n] = 3e0 * c[n]
      d[n] = d[n-1]
      return
end
