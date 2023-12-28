
using .cppdefs


function findzall()
  local dztop = @fortGet("dztop", rc_kind)
  local dzbot = cppdefs.fixed_bottom_thickness ? @fortGet("dzbot", rc_kind) : 9090909090e0
  local dnkm1 = rc_kind(NK-1)
  local dnkm1inv = 1e0 / dnkm1
  
  # stretching of the vertical grid
  #    pfac is the stretching in z. higher pfac gives more points near surf.
  local pfac = 2e0
  
  local epm1 = exp(pfac) - 1e0
  local epm1inv = 1e0 / epm1
  
  for j in 0:NJ+1
    for i in 0:NI+1
      local hpd = h[i+1, j+1] * @fortGet("hdl",rc_kind) + @fortGet("dztop", rc_kind)
      for k in NK:NK+1
        local sigma = rc_kind(k) - 0.5e0
        zc[i, j, k] = (sigma - dnkm1) * hpd - @fortGet("dztop", rc_kind)
      end
      for k in NK-1:NK+1
        local sigma = rc_kind(k)
        zf[i, j, k] = (sigma - dnkm1) * hpd - @fortGet("dztop", rc_kind)
      end
    end
  end
  
  
  @static if (cppdefs.fixed_bottom_thickness)
    dnkm1 = rc_kind(NK - 2)
    dnkm1inv = 1e0 / dnkm1
  end
  
  for j in 0:NJ+1
    for i in 0:NI+1
      for k in 0:NK-1
        @static if (cppdefs.fixed_bottom_thickness)
          local sigma = rc_kind(k-1) - 0.5e0
        else
          local sigma = rc_kind(k) - 0.5e0
        end
        local xfac = (dnkm1 - sigma) * dnkm1inv
        
        @static if (cppdefs.fixed_bottom_thickness)
          zc[i, j, k] = (exp(pfac*xfac) - 1e0) * epm1inv * (D[i,j] + dzbot + dztop) - dztop
        else
          zc[i, j, k] = (exp(pfac*xfac) - 1e0) * epm1inv * (D[i,j] + dztop) - dztop
        end

      end
      for k in -1:NK-2
        @static if (cppdefs.fixed_bottom_thickness)
          local sigma = rc_kind(k-1)
        else
          local sigma = rc_kind(k)
        end
        local xfac = (dnkm1 - sigma) * dnkm1inv
        @static if (cppdefs.fixed_bottom_thickness)
          zf[i, j, k] = (exp(pfac*xfac) - 1e0) * epm1inv * (D[i,j] + dzbot + dztop) - dztop
        else
          zf[i, j, k] = (exp(pfac*xfac) - 1e0) * epm1inv * (D[i,j] + dztop) - dztop
        end
      end
      
      @static if (cppdefs.fixed_bottom_thickness)
        for k in 0:1
          local sigma = rc_kind(k) - 0.5e0
          zc[i, j, k] = ((sigma - 0e0) / 1e0) * dzbot + D[i, j]
        end
        for k in -1:1
          local sigma = rc_kind(k)
          zf[i, j, k] = ((sigma - 0e0) / 1e0) * dzbot + D[i, j]
        end
      end
      
    end
  end

end