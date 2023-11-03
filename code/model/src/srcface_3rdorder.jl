

function srcface(n, step)
  local uxi = zeros(rc_kind, (NI+1, NJ))
  local uyi = zeros(rc_kind, (NI+1, NJ))
  local vxj = zeros(rc_kind, (NI, NJ+1))
  local vyj = zeros(rc_kind, (NI, NJ+1))


  local fac = EPS * @fortGet("delta", rc_kind)
  local ainv = 1e0 / apr
  local fac2 = EPS * @fortGet("lambda", rc_kind)
  local be2 = @fortGet("beta", rc_kind) * EPS * EPS

  for j in 1:NJ
    for i in 1:NI-1
      uxi[i+1, j] = 0.5e0 * (ux[i+2, j+1] + ux[i+1, j+1])
      uyi[i+1, j] = 0.5e0 * (uy[i+2, j+1] + uy[i+1, j+1])
    end
    uxi[NI+1, j] = 0.5e0 * (ux[2, j+1] + ux[NI+1, j+1])
    uyi[NI+1, j] = 0.5e0 * (uy[2, j+1] + uy[NI+1, j+1])
    uxi[1, j] = (uxi[NI+1, j])
    uyi[1, j] = (uyi[NI+1, j])
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
        local uint = xa * (u[i0+1, j+1, k+1, n+1] + u[ip1+1, j+1, k+1, n+1]) + xb * (u[im1+1, j+1, k+1, n+1] + u[ip2+1, j+1, k+1, n+1])
        local vint = xa * (v[i0+1, j+1, k+1, n+1] + v[ip1+1, j+1, k+1, n+1]) + xb * (v[im1+1, j+1, k+1, n+1] + v[ip2+1, j+1, k+1, n+1])
        local wint = xa * (w[i0+1, j+1, k+1, n+1] + w[ip1+1, j+1, k+1, n+1]) + xb * (w[im1+1, j+1, k+1, n+1] + w[ip2+1, j+1, k+1, n+1])
        
        local vcif = EPS * (xa * (uvis[i0, j, k] + uvis[ip1, j, k]) + xb * (uvis[im1, j, k] + uvis[ip2, j, k]))
        local vcjf = EPS * (xa * (vvis[i0, j, k] + vvis[ip1, j, k]) + xb * (vvis[im1, j, k] + vvis[ip2, j, k])) 
        
        
        local px = ((p[ip1+1, j+1, k+1] - p[i+1, j+1, k+1]) * gqi[i+1, j, k, 1] 
        + 0.25e0 * (p[ip1+1, j+2, k+1] + p[i+1, j+2, k+1] - p[ip1+1, j, k+1] - p[i+1, j, k+1]) * gqi[i+1, j, k, 2]
        + 0.25e0 * (p[ip1+1, j+1, k+2] + p[i+1, j+1, k+2] - p[ip1+1, j+1, k] - p[i+1, j+1, k]) * gqi3[i+1, j, k])
        sifc[i+1, j, k] = ((uxi[i+1, j] * (-ffi[i+1, j] * vint + fac * bbi[i+1, j] * wint - vcif) 
        + uyi[i+1, j] * ( ffi[i+1, j] * uint - vcjf)) * Jifc[i+1, j, k] + grpifc[i+1, j, k]) + px
        
      end
    end
  end
  
  for k in 1:NK
    for j in 1:NJ
      sifc[1,j,k] = sifc[NI+1, j ,k]
    end
  end
  
  for j in 0:NJ
    for i in 1:NI
      vxj[i, j+1] = 0.5e0 * (vx[i+1, j+2] + vx[i+1, j+1])
      vyj[i, j+1] = 0.5e0 * (vy[i+1, j+2] + vy[i+1, j+1])
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
        local vint = xa * (v[i+1, j0+1, k+1, n+1] + v[i+1, jp1+1, k+1, n+1]) + xb * (v[i+1, jm1+1, k+1, n+1] + v[i+1, jp2+1, k+1, n+1])
        local uint = xa * (u[i+1, j0+1, k+1, n+1] + u[i+1, jp1+1, k+1, n+1]) + xb * (u[i+1, jm1+1, k+1, n+1] + u[i+1, jp2+1, k+1, n+1])
        local wint = xa * (w[i+1, j0+1, k+1, n+1] + w[i+1, jp1+1, k+1, n+1]) + xb * (w[i+1, jm1+1, k+1, n+1] + w[i+1, jp2+1, k+1, n+1])
        local vcif = EPS * (xa * (uvis[i, j0, k] + uvis[i, jp1, k]) + xb * (uvis[i, jm1, k] + uvis[i, jp2, k]))
        local vcjf = EPS * (xa * (vvis[i, j0, k] + vvis[i, jp1, k]) + xb * (vvis[i, jm1, k] + vvis[i, jp2, k]))
        
        local py = ((p[i+1, j+2, k+1] - p[i+1, j+1, k+1]) * gqj[i, j+1, k, 2]
        + 0.25e0 * (p[i+2, j+2, k+1] + p[i+2, j+1, k+1] - p[i, j+2, k+1] - p[i, j+1, k+1]) * gqj[i, j+1, k, 1]
        + 0.25e0 * (p[i+1, j+2, k+2] + p[i+1, j+1, k+2] - p[i+1, j+2, k] - p[i+1, j+1, k]) * gqj3[i, j+1, k])
        
        sjfc[i, j+1, k] = ((vxj[i, j+1] * (-ffj[i, j+1] * vint + fac * bbj[i, j+1] * wint - vcif)
        + vyj[i, j+1] * ( ffc[i+1, j+1] * uint - vcjf)) * Jjfc[i, j+1, k] + grpjfc[i, j+1, k]) + py
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

        local vint = 0.5e0 * (3e0 * v[i+1, in+1, k+1, n+1] - v[i+1, inm1+1, k+1, n+1])
        local uint = 0.5e0 * (3e0 * u[i+1, in+1, k+1, n+1] - u[i+1, inm1+1, k+1, n+1])
        local wint = 0.5e0 * (3e0 * w[i+1, in+1, k+1, n+1] - w[i+1, inm1+1, k+1, n+1])

        local vcif= EPS*uvis[i,in,k] 
        local vcjf= EPS*vvis[i,in,k]
        
        local py = ((p[i+1, j+2, k+1] - p[i+1, j+1, k+1]) * gqj[i, j+1, k, 2]
          + 0.25e0 * (p[i+2, j+2, k+1] + p[i+2, j+1, k+1] - p[i, j+2, k+1] - p[i, j+1, k+1]) * gqj[i, j+1, k, 1]
          + 0.25e0 * (p[i+1, j+2, k+2] + p[i+1, j+1, k+2] - p[i+1, j+2, k] - p[i+1, j+1, k]) * gqj3[i, j+1, k])

        sjfc[i, j+1, k] = ((vxj[i, j+1] * (-ffj[i, j+1] * vint + fac * bbj[i, j+1] * wint - vcif)
        + vyj[i, j+1] * ( ffc[i+1, j+1] * uint - vcjf)) * Jjfc[i, j+1, k] + grpjfc[i, j+1, k]) + py
      end
    end
  end


end

