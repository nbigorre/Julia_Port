function ini_setup(pcorr)

  #     subroutine init(p,vfent,ufex)                                    
  #     ------------------------------------------------                  
  #     For the curvilinear grid.                                         
  #     sigma levels are evenly spaced.  Hence wz is a function of        
  #     time but not of z.   Here u,v,w refer to the xi,eta,sigma direcito
  #     Those metric quantities that do not change with time are evaluated
  #     The rest are evaluated in "sigma.f" which will be called at every 
  #     step.                                                             
  #     Also it is absolutely essential that the integral of the flux     
  #     over the entrance section is equal to that over the exit section. 
  #     -No longer necessary with the free surface.                       
  #     This may require some special attention when the free surface is  
  #     moving at the exit. At the entrance vf is held fixed.             
  #                                                                       

  #REAL(kind=rc_kind) ::   xdu(0:NI+1,0:NJ+1),ydu(0:NI+1,0:NJ+1),xdv(0:NI+1,0:NJ+1),ydv(0:NI+1,0:NJ+1)
  #REAL(kind=rc_kind) :: pcorr(maxout)
  #REAL(kind=rc_kind) :: fconst,phi0,cosphi0,sinphi0,dthet,dtheta,dphi,dep,y
  #INTEGER i,j,k,jmid,iseed 
  local xdu = OffsetArrays.zeros(rc_kind, 0:NI+1, 0:NJ+1)
  local ydu = OffsetArrays.zeros(rc_kind, 0:NI+1, 0:NJ+1)
  local xdv = OffsetArrays.zeros(rc_kind, 0:NI+1, 0:NJ+1)
  local ydv = OffsetArrays.zeros(rc_kind, 0:NI+1, 0:NJ+1)
  local DLinv = 1e0 / DL
  local phi0 = phi0deg * PI / 180e0  #     phi0 is the central lat, dphi,dtheta grid spacing(angle)          
  local dphi = dy / (apr * AL)



  #--------------------------------------------------
  # INITIALIZATION OF THE GRID
  #


  # HORIZONTAL GRID

  xc[0] = -0.5 * dx * 1e-3
  yc[0] = -0.5 * dy * 1e-3

  for i in 1:NI+1
    xc[i] = xc[i-1] + dx * 1e-3
  end
  for j in 1:NJ+1
    yc[j] = yc[j-1] + dy * 1e-3
  end

  for j in 0:NJ+1
    for i in 0:NI+1
      xdu[i, j] = dx / LEN
      ydu[i, j] = 0e0
      xdv[i, j] = 0e0
      ydv[i, j] = dy / LEN # used for constant dy
      # ydv[i,j]= dyM[j]/LEN 
    end
  end

  for j in 0:NJ+1
    for i in 0:NI+1
      J2d[i, j] = xdu[i, j] * ydv[i, j] - xdv[i, j] * ydu[i, j]
      ux[i, j] = ydv[i, j] / J2d[i, j]
      vx[i, j] = -ydu[i, j] / J2d[i, j]
      uy[i, j] = -xdv[i, j] / J2d[i, j]
      vy[i, j] = xdu[i, j] / J2d[i, j]
      g11[i, j] = ux[i, j] * ux[i, j] + uy[i, j] * uy[i, j]
      g12[i, j] = ux[i, j] * vx[i, j] + uy[i, j] * vy[i, j]
      g22[i, j] = vx[i, j] * vx[i, j] + vy[i, j] * vy[i, j]
    end
  end


  # DEPTH OF EVERY MODEL COLUMN

  if (lv_flat_bottom == 0)
    ini_topog()
    #@ccall "./PSOM_LIB.so".ini_topog_()::Cvoid
  else
    local dep = total_depth
    @views @. D[:, :] .= -dep * DLinv                             # D(i,j) is -ve and  non-dim by DL
    smooth()
    #@ccall "./PSOM_LIB.so".smooth_()::Cvoid          # Compute partial derivatives dD/dx,dD/dy 
    # Does not smooth anything.   
  end

  #----------------------------------------------
  # COMPUTATION OF THE PLANETARY VORTICITY.
  #


  local fconst = 2e0 * OMEGA / FPAR
  local jmid = NJ / 2
  for j in 0:NJ+1
    latrad[j] = phi0 + rc_kind(1 - fplane) * (Float64(j - jmid) * dphi)
    #latrad[j] = phi0
  end

  for j in 1:NJ
    @views @. ffi[:, j] .= fconst * sin(latrad[j])
    @views @. bbi[:, j] .= fnhhy * fconst * cos(latrad[j])
  end
  for j in 0:NJ
    @views @. ffj[:, j] .= fconst * sin(0.5e0 * (latrad[j+1] + latrad[j]))
    @views @. bbj[:, j] .= fnhhy * fconst * cos(0.5e0 * (latrad[j+1] + latrad[j]))
  end
  for j in 0:NJ+1
    @views @. ffc[:, j] .= fconst * sin(latrad[j])
    @views @. bbc[:, j] .= fnhhy * fconst * cos(latrad[j])
  end



  #--------------------------------------------------
  # INITIALIZATION OF THE VARIABLES
  #

  @views @. u[:, :, :, 0] .= 0e0
  @views @. v[:, :, :, 0] .= 0e0
  @views @. w[:, :, :, 0] .= 0e0
  @views @. s[:, :, :, 0] .= 0e0
  @views @. T[:, :, :, 0] .= 0e0
  @views @. Tr[:, :, :, :, 0] .= 0e0
  @views @. conv[:, :, :] .= 0
  @views @. con100[:, :, :] .= 0

  @views @. pcorr[:] .= 0e0

  @views @. uvis[:, :, :] .= 0e0
  @views @. vvis[:, :, :] .= 0e0
  @views @. wvis[:, :, :] .= 0e0

  h .= 0e0
  uf .= 0e0
  vf .= 0e0
  wf .= 0e0

  @views @. ufbce[1:NJ, 1:NK] = uf[NI,1:NJ, 1:NK]
  @views @. ufbcw[1:NJ, 1:NK] = uf[0,1:NJ, 1:NK]


  @views @. vfbcn[:, :] .= vf[:, NJ, :]
  @views @. vfbcs[:, :] .= vf[:, 0, :]

  @views @. wfbcb[:, :] .= 0e0    # at the sea bed, wf=0                                              

  #---------------------

  findzall()
  #@ccall "./PSOM_LIB.so".findzall_()::Cvoid                                            # finds the vertical grid

  ini_st()
  #@ccall "./PSOM_LIB.so".ini_st_()::Cvoid                                              # initialize s,T                                                    
  
  evalrho(rho, 0)
  #@ccall "./PSOM_LIB.so".evalrho_(pointer(rho)::Ptr{rc_kind}, Ref(0)::Ref{Int})::Cvoid       # deduce rho from s,T 

  ini_h()
  #@ccall "./PSOM_LIB.so".ini_h_()::Cvoid                                               # initialize free surface h

  findzall()
  #@ccall "./PSOM_LIB.so".findzall_()::Cvoid                                            # find z again


  @static if (cppdefs.file_output)
    if (lv_flat_bottom != 0)

      local fs = open(string(dirout, "/zgrid.out"), "w")
      write(fs, "#vertical grid\n")
      for k in 0:NK+1
        write(fs, string(k, "  ", zc[10, 10, k] * 1000e0, "\n"))
      end
      write(fs, "# face values\n")
      for k in -1:NK+1
        write(fs, string(k, "  ", zf[10, 10, k] * 1000e0, "\n"))
      end
    end

  end
  #endif

  @views @. advecpv[:] .= 0e0
  @views @. friction[:] .= 0e0
  @views @. diabatic[:] .= 0e0

end
