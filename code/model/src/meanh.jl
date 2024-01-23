function meanh(NI, NJ, h)
  #     -------------------------------                                   
  #     computes the mean free-surface elevation within the domain        
  # use header, only : rc_kind
  # integer NI,NJ,i,j 
  # REAL(kind=rc_kind) :: h(0:NI+1,0:NJ+1),hsum,hmean 

  local hsum = 0e0
  for i in 1:NI
    for j in 1:NJ
      hsum = hsum + h[i, j]
    end
  end

  local hmean = hsum / rc_kind(NI * NJ)

  if (abs(hmean) > 100e0)
    println("error, hmean=", hmean)
    exit(1)
  end

  return hmean

end