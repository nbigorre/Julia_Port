function seval(n, u, x, y, b, c, d)

  #-----------------------------------------------------------            
  #     use header, only : rc_kind
  #  integer n 
  #  REAL(kind=rc_kind) :: seval
  #  REAL(kind=rc_kind) ::  u, x(n), y(n), b(n), c(n), d(n) 
  #                                                                       
  #  this subroutine evaluates the cubic spline function                  
  #                                                                       
  #    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3 
  #                                                                       
  #    where  x(i) .lt. u .lt. x(i+1), using horner's rule                
  #                                                                       
  #  if  u .lt. x(1) then  i = 1  is used.                                
  #  if  u .ge. x(n) then  i = n  is used.                                
  #                                                                       
  #  input..                                                              
  #                                                                       
  #    n = the number of data points                                      
  #    u = the abscissa at which the spline is to be evaluated            
  #    x,y = the arrays of data abscissas and ordinates                   
  #    b,c,d = arrays of spline coefficients computed by spline           
  #                                                                       
  #  if  u  is not in the same interval as the previous call, then a      
  #  binary search is performed to determine the proper interval.         
  #                                                                       
  #  integer i, j, k 
  #  REAL(kind=rc_kind) :: dx 
  local i = 1
  if (i >= n)
    i = 1
  end
  if ((u < x[i]))
    @goto g10
  end
  if (u <= x[i+1])
    @goto g30
  end
  #                                                                       
  #  binary search                                                        
  #                                                                       
  @label g10
  i = 1
  local j = n + 1

  @label g20
  local k = div(i + j, 2)
  if (u < x[k])
    j = k
  end
  if (u >= x[k])
    i = k
  end
  if (j > i + 1)
    @goto g20
  end
  #                                                                       
  #  evaluate spline                                                      
  #   
  @label g30                                                                    
  local dx = u - x[i]
  return y[i] + dx * (b[i] + dx * (c[i] + dx * d[i]))
end
