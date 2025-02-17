

function srcface(n, step)
  local uxi = OffsetArray(zeros(rc_kind, (NI+1, NJ)), 0:NI, 1:NJ)
  local uyi = OffsetArray(zeros(rc_kind, (NI+1, NJ)), 0:NI, 1:NJ)
  local vxj = OffsetArray(zeros(rc_kind, (NI, NJ+1)), 1:NI, 0:NJ)
  local vyj = OffsetArray(zeros(rc_kind, (NI, NJ+1)), 1:NI, 0:NJ)


  local fac = EPS * delta
  local ainv = 1e0 / apr
  local fac2 = EPS * lambda
  local be2 = beta * EPS * EPS

  for j in 1:NJ
    for i in 1:NI-1
      uxi[i, j] = 0.5e0 * (ux[i+1, j] + ux[i, j])
      uyi[i, j] = 0.5e0 * (uy[i+1, j] + uy[i, j])
    end
    uxi[NI, j] = 0.5e0 * (ux[1, j] + ux[NI, j])
    uyi[NI, j] = 0.5e0 * (uy[1, j] + uy[NI, j])
    uxi[0, j] = (uxi[NI, j])
    uyi[0, j] = (uyi[NI, j])
  end
  
  for i in 1:NI
    local xa::rc_kind = 9e0 / 16e0
    local xb::rc_kind = -1e0 / 16e0
    
    local i0::Int = i 
    local im1::Int = i-1 
    local ip1::Int = i+1 
    local ip2::Int = i+2
    
    if (i == 1)
      i0 = i 
      im1 = NI 
      ip1 = i+1 
      ip2 = i+2
    elseif (i == NI - 1)
      i0 = i 
      im1 = i-1 
      ip1 = i+1 
      ip2 = 1 
    elseif (i == NI) 
      i0 = i 
      im1 = i-1 
      ip1 = 1 
      ip2 = 2 
    end
    
    for k in 1:NK
      for j in 1:NJ
        local uint = xa * (u[i0, j, k, n] + u[ip1, j, k, n]) + xb * (u[im1, j, k, n] + u[ip2, j, k, n])
        local vint = xa * (v[i0, j, k, n] + v[ip1, j, k, n]) + xb * (v[im1, j, k, n] + v[ip2, j, k, n])
        local wint = xa * (w[i0, j, k, n] + w[ip1, j, k, n]) + xb * (w[im1, j, k, n] + w[ip2, j, k, n])
        
        local vcif = EPS * (xa * (uvis[i0, j, k] + uvis[ip1, j, k]) + xb * (uvis[im1, j, k] + uvis[ip2, j, k]))
        local vcjf = EPS * (xa * (vvis[i0, j, k] + vvis[ip1, j, k]) + xb * (vvis[im1, j, k] + vvis[ip2, j, k])) 
        
        
        local px = ((p[ip1, j, k] - p[i, j, k]) * gqi[i, j, k, 1] 
        + 0.25e0 * (p[ip1, j+1, k] + p[i, j+1, k] - p[ip1, j-1, k] - p[i, j-1, k]) * gqi[i, j, k, 2]
        + 0.25e0 * (p[ip1, j, k+1] + p[i, j, k+1] - p[ip1, j, k-1] - p[i, j, k-1]) * gqi3[i, j, k])
        sifc[i, j, k] = ((uxi[i, j] * (-ffi[i, j] * vint + fac * bbi[i, j] * wint - vcif) 
        + uyi[i, j] * ( ffi[i, j] * uint - vcjf)) * Jifc[i, j, k] + grpifc[i, j, k]) + px
        
      end
    end
  end
  
  for k in 1:NK
    for j in 1:NJ
      sifc[0,j,k] = sifc[NI, j ,k]
    end
  end
  
  for j in 0:NJ
    for i in 1:NI
      vxj[i, j] = 0.5e0 * (vx[i, j+1] + vx[i, j])
      vyj[i, j] = 0.5e0 * (vy[i, j+1] + vy[i, j])
    end
  end
  
  
  for j in 1:NJ-1
    local j0= j 
    local jm1= j-1 
    local jp1= j+1 
    local jp2= j+2 
    local xa= 9e0/16e0 
    local xb= -1e0/16e0
    
    if (j == 1)
      j0= j 
      jm1=NJ 
      jp1= j+1 
      jp2= j+2 
      xa = 0.5e0 
      xb = 0e0 
    elseif (j == NJ-1)
      j0= j 
      jm1= j-1 
      jp1= NJ 
      jp2= 1 
      xa = 0.5e0 
      xb = 0e0
    end
    
    for k in 1:NK
      for i in 1:NI
        local vint = xa * (v[i, j0, k, n] + v[i, jp1, k, n]) + xb * (v[i, jm1, k, n] + v[i, jp2, k, n])
        local uint = xa * (u[i, j0, k, n] + u[i, jp1, k, n]) + xb * (u[i, jm1, k, n] + u[i, jp2, k, n])
        local wint = xa * (w[i, j0, k, n] + w[i, jp1, k, n]) + xb * (w[i, jm1, k, n] + w[i, jp2, k, n])
        local vcif = EPS * (xa * (uvis[i, j0, k] + uvis[i, jp1, k]) + xb * (uvis[i, jm1, k] + uvis[i, jp2, k]))
        local vcjf = EPS * (xa * (vvis[i, j0, k] + vvis[i, jp1, k]) + xb * (vvis[i, jm1, k] + vvis[i, jp2, k]))
        
        local py = ((p[i, j+1, k] - p[i, j, k]) * gqj[i, j, k, 2]
        + 0.25e0 * (p[i+1, j+1, k] + p[i+1, j, k] - p[i-1, j+1, k] - p[i-1, j, k]) * gqj[i, j, k, 1]
        + 0.25e0 * (p[i, j+1, k+1] + p[i, j, k+1] - p[i, j+1, k-1] - p[i, j, k-1]) * gqj3[i, j, k])
        
        sjfc[i, j, k] = ((vxj[i, j] * (-ffj[i, j] * vint + fac * bbj[i, j+1] * wint - vcif)
        + vyj[i, j] * ( ffc[i, j] * uint - vcjf)) * Jjfc[i, j, k] + grpjfc[i, j, k]) + py
      end
    end
  end
  
  for k in 1:NK
    for i in 1:NI
      for j in 0:NJ:NJ
        local in=1 
        local inm1= 2 
        if (j == NJ)
          in = NJ
          inm1 = NJ-1
        end

        local vint = 0.5e0 * (3e0 * v[i, in, k, n] - v[i, inm1, k, n])
        local uint = 0.5e0 * (3e0 * u[i, in, k, n] - u[i, inm1, k, n])
        local wint = 0.5e0 * (3e0 * w[i, in, k, n] - w[i, inm1, k, n])

        local vcif= EPS*uvis[i,in,k] 
        local vcjf= EPS*vvis[i,in,k]
        
        local py = ((p[i, j+1, k] - p[i, j, k]) * gqj[i, j, k, 2]
          + 0.25e0 * (p[i+1, j+1, k] + p[i+1, j, k] - p[i-1, j+1, k] - p[i-1, j, k]) * gqj[i, j, k, 1]
          + 0.25e0 * (p[i, j+1, k+1] + p[i, j, k+1] - p[i, j+1, k-1] - p[i, j, k-1]) * gqj3[i, j, k])

        sjfc[i, j, k] = ((vxj[i, j] * (-ffj[i, j] * vint + fac * bbj[i, j] * wint - vcif)
        + vyj[i, j] * ( ffc[i, j] * uint - vcjf)) * Jjfc[i, j, k] + grpjfc[i, j, k]) + py
      end
    end
  end


end

