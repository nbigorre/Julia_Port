using OffsetArrays


global const rpevalgrad_Sche_dZ = OffsetArray(zeros(rc_kind, (NI + 2, NK + 1)), 0:NI+1, 0:NK)
global const rpevalgrad_Sche_dR = OffsetArray(zeros(rc_kind, (NI + 2, NK + 1)), 0:NI+1, 0:NK)
global const rpevalgrad_Sche_FC = OffsetArray(zeros(rc_kind, (NI + 2, NJ + 2)), 0:NI+1, 0:NJ+1)
global const rpevalgrad_Sche_dZx = OffsetArray(zeros(rc_kind, (NI + 2, NJ + 2)), 0:NI+1, 0:NJ+1)
global const rpevalgrad_Sche_dRx = OffsetArray(zeros(rc_kind, (NI + 2, NJ + 2)), 0:NI+1, 0:NJ+1)
global const rpevalgrad_Sche_rx = OffsetArray(zeros(rc_kind, (NI + 2, NJ + 2)), 0:NI+1, 0:NJ+1)
global const rpevalgrad_Sche_dn_u = OffsetArray(zeros(rc_kind, (NI + 2, NJ + 2)), 0:NI+1, 0:NJ+1)
global const rpevalgrad_Sche_dm_v = OffsetArray(zeros(rc_kind, (NI + 2, NJ + 2)), 0:NI+1, 0:NJ+1)
global const rpevalgrad_Sche_rhol = OffsetArray(zeros(rc_kind, (NI + 2, NJ + 2, NK + 2)), 0:NI+1, 0:NJ+1, 0:NK+1)
global const rpevalgrad_Sche_z_r = OffsetArray(zeros(rc_kind, (NI + 2, NJ + 2, NK + 2)), 0:NI+1, 0:NJ+1, 0:NK+1)
global const rpevalgrad_Sche_Hz = OffsetArray(zeros(rc_kind, (NI + 2, NJ + 2, NK + 2)), 0:NI+1, 0:NJ+1, 0:NK+1)
global const rpevalgrad_Sche_z_w = OffsetArray(zeros(rc_kind, (NI + 2, NJ + 2, NK + 3)), 0:NI+1, 0:NJ+1, -1:NK+1)
global const rpevalgrad_Sche_P_Sche = OffsetArray(zeros(rc_kind, (NI + 2, NJ + 2, NK)), 0:NI+1, 0:NJ+1, 1:NK)



function rpevalgrad_Sche(nl::Int)

  local dZ = rpevalgrad_Sche_dZ 
  local dR = rpevalgrad_Sche_dR 
  local FC = rpevalgrad_Sche_FC 
  local dZx = rpevalgrad_Sche_dZx 
  local dRx = rpevalgrad_Sche_dRx 
  local rx = rpevalgrad_Sche_rx 
  local dn_u = rpevalgrad_Sche_dn_u 
  local dm_v = rpevalgrad_Sche_dm_v 
  local rhol = rpevalgrad_Sche_rhol 
  local z_r = rpevalgrad_Sche_z_r 
  local Hz = rpevalgrad_Sche_Hz 
  local z_w = rpevalgrad_Sche_z_w 
  local P_Sche = rpevalgrad_Sche_P_Sche 

  local epsil = 0e0


  # ====================================================================================h
  #
  # -----------------------------------------------------------------
  # Part I: The variables of PSOM are translated into ROMS variables.
  # -----------------------------------------------------------------
  #
  #
  # --------
  # Indices

  local imin = 1
  local imax = NI

  local jmin = 1
  local jmax = NJ

  local N = NK

  local istrU = imin
  local iend = imax
  local jstrV = jmin
  local jend = jmax
  local istr = imin
  local jstr = jmin

  # --------
  # Geometry

  dn_u .= dx
  dm_v .= dy

  @views @. z_r = DL * zc
  @views @. z_w = DL * zf

  for k in 1:NK
    Hz[0:NI+1, 0:NJ+1, k] = (z_w[:, :, k] .- z_w[:, :, k-1])
  end

  # ----------
  # Geophysics

  local g = gpr * 10e0
  local rho0 = R0

  evalrho(rho, nl)
  rhol .= rho .- rho0

  # ------------------------
  # Boundary condition in x:

  Hz[0, :, :] = Hz[NI, :, :]
  Hz[NI+1, :, :] = Hz[1, :, :]

  rhol[0, :, :] = rhol[NI, :, :]
  rhol[NI+1, :, :] = rhol[1, :, :]



  # ====================================================================================h

  # -----------------------------------------------------------------
  # Part II: The code.
  # -----------------------------------------------------------------
  # The cpp tests have been removes.
  # Except from explicitely designed ligns, nothing has been changed.

  #
  # Preliminary step (same for XI- and ETA-components):
  #------------ ---- ----- --- --- --- ----------------
  #

  local GRho = g / rho0
  local HalfGRho = 0.5e0 * GRho

  for j in jstrV-1:jend
    for k in 1:N-1
      for i in istrU-1:iend
        dZ[i, k] = z_r[i, j, k+1] - z_r[i, j, k]
        dR[i, k] = rhol[i, j, k+1] - rhol[i, j, k]
      end
    end
    for i in istrU-1:iend
      dR[i, N] = dR[i, N-1]
      dR[i, 0] = dR[i, 1]
      dZ[i, N] = dZ[i, N-1]
      dZ[i, 0] = dZ[i, 1]
    end
    for k in N:-1:1
      for i in istrU-1:iend
        local cff = 2e0 * dZ[i, k] * dZ[i, k-1]
        dZ[i, k] = cff / (dZ[i, k] + dZ[i, k-1])

        local cfr = 2e0 * dR[i, k] * dR[i, k-1]
        if (cfr > epsil)
          dR[i, k] = cfr / (dR[i, k] + dR[i, k-1])
        else
          dR[i, k] = 0e0
        end
      end
    end
    for i in istrU-1:iend
      P_Sche[i, j, N] = 0e0 * g * z_w[i, j, N] + GRho * (rhol[i, j, N] + 0.5 * (rhol[i, j, N] - rhol[i, j, N-1]) * (z_w[i, j, N] - z_r[i, j, N]) / (z_r[i, j, N] - z_r[i, j, N-1])) * (z_w[i, j, N] - z_r[i, j, N])
    end
    for k in N-1:-1:1
      for i in istrU-1:iend
        P_Sche[i, j, k] = P_Sche[i, j, k+1] + HalfGRho * ((rhol[i, j, k+1] + rhol[i, j, k]) * (z_r[i, j, k+1] - z_r[i, j, k])
                                                                  -
                                                                  0.2e0 * ((dR[i, k+1] - dR[i, k]) * (z_r[i, j, k+1] - z_r[i, j, k] - (1e0 / 12e0) * (dZ[i, k+1] + dZ[i, k]))
                                                                           -
                                                                           (dZ[i, k+1] - dZ[i, k]) * (rhol[i, j, k+1] - rhol[i, j, k] - (1e0 / 12e0) * (dR[i, k+1] + dR[i, k]))))
      end
    end
  end


  # Caution:
  # The line P_Sche(i,j,N)=... has been modified !
  # The term g*z_w(i,j,N) has been set to 0 because PSOM only requires the baroclinic terms and not 
  #  the barotropic part linked with surface elevation (g.R0.h).
  #


  #
  # Compute XI-component of pressure gradient term:
  #-------- ------------ -- -------- -------- -----
  #
  for k in N:-1:1
    for j in jstr:jend
      for i in imin:imax
        FC[i, j] = z_r[i, j, k] - z_r[i-1, j, k]
        rx[i, j] = rhol[i, j, k] - rhol[i-1, j, k]
      end
    end

    # This part has been included to take into account the boundary conditions:
    #++++ Added by JBG, 20120425.
    @. @views FC[0, :] = FC[NI, :]
    @. @views FC[1, :] = z_r[1, :, k] - z_r[NI, :, k]
    @. @views rx[0, :] = rx[NI-1, :]
    @. @views rx[1, :] = rhol[1, :, k] - rhol[NI, :, k]
    @. @views FC[NI+1, :] = FC[1, :]
    @. @views rx[NI+1, :] = rx[1, :]

    for j in jstr:jend
      for i in istrU-1:iend
        local cff = 2e0 * FC[i, j] * FC[i+1, j]
        if (cff > epsil)
          dZx[i, j] = cff / (FC[i, j] + FC[i+1, j])
        else
          dZx[i, j] = 0e0
        end

        local cfr = 2e0 * rx[i, j] * rx[i+1, j]
        if (cfr > epsil)
          dRx[i, j] = cfr / (rx[i, j] + rx[i+1, j])
        else
          dRx[i, j] = 0e0
        end
      end

      for i in istrU:iend
        ru_Sche[i, j, k] = 0.5e0 * (Hz[i, j, k] + Hz[i-1, j, k]) * dn_u[i, j] * (
                                 (P_Sche[i-1, j, k] - P_Sche[i, j, k]) - HalfGRho * (
                                   (rhol[i, j, k] + rhol[i-1, j, k]) * (z_r[i, j, k] - z_r[i-1, j, k])
                                   -
                                   0.2e0 * ((dRx[i, j] - dRx[i-1, j]) * (z_r[i, j, k] - z_r[i-1, j, k]
                                                                             -
                                                                             (1e0 / 12e0) * (dZx[i, j] + dZx[i-1, j]))
                                            -
                                            (dZx[i, j] - dZx[i-1, j]) * (rhol[i, j, k] - rhol[i-1, j, k]
                                                                             -
                                                                             (1e0 / 12e0) * (dRx[i, j] + dRx[i-1, j]))
                                   )
                                 )
                               )
      end
    end

    #
    # ETA-component of pressure gradient term:
    #-------------- -- -------- -------- -----
    #


    for j in jmin:jend
      for i in istr:iend
        FC[i, j] = (z_r[i, j, k] - z_r[i, j-1, k])
        rx[i, j] = (rhol[i, j, k] - rhol[i, j-1, k])
      end
    end

    # This part has been altered to take into account the boundary conditions:
    #++++ Modified by JBG, 20120425

    for i in istr:iend
      FC[i, jmin] = FC[i, jmin+1]
      FC[i, jmin-1] = FC[i, jmin]
      rx[i, jmin] = rx[i, jmin+1]
      rx[i, jmin-1] = rx[i, jmin]
    end

    for i in istr:iend
      FC[i, jmax] = FC[i, jmax-1]
      FC[i, jmax+1] = FC[i, jmax]
      rx[i, jmax] = rx[i, jmax-1]
      rx[i, jmax+1] = rx[i, jmax-1]
    end


    for j in jstrV-1:jend
      for i in istr:iend
        local cff = 2e0 * FC[i, j] * FC[i, j+1]
        if (cff > epsil)
          dZx[i, j] = cff / (FC[i, j] + FC[i, j+1])
        else
          dZx[i, j] = 0e0
        end

        local cfr = 2e0 * rx[i, j] * rx[i, j+1]
        if (cfr > epsil)
          dRx[i, j] = cfr / (rx[i, j] + rx[i, j+1])
        else
          dRx[i, j] = 0e0
        end
      end

      if (j >= jstrV)
        for i in istr:iend
          rv_Sche[i, j, k] = 0.5e0 * (Hz[i, j, k] + Hz[i, j-1, k]) * dm_v[i, j] * (
                                   (P_Sche[i, j-1, k] - P_Sche[i, j, k]) - HalfGRho * (
                                     (rhol[i, j, k] + rhol[i, j-1, k]) * (z_r[i, j, k] - z_r[i, j-1, k]) - 0.2e0 * (
                                       (dRx[i, j] - dRx[i, j-1]) * (z_r[i, j, k] - z_r[i, j-1, k] - (1e0 / 12e0) * (dZx[i, j] + dZx[i, j-1]))
                                       -
                                       (dZx[i, j] - dZx[i, j-1]) * (rhol[i, j, k] - rhol[i, j-1, k] - (1e0 / 12e0) * (dRx[i, j] + dRx[i, j-1]))
                                     )
                                   )
                                 )
        end
      end
    end

  end



  # ====================================================================================h

  # -----------------------------------------------------------------
  # Part III: The variables of ROMS are now translated into PSOM variables.
  # -----------------------------------------------------------------

  #--------------------------------
  # ru_Sche and rv_Sche are the equivalent to grpifc and grpjfc but 
  # - are dimensionalized
  # - have a shift of one grid point.
  # - the boundary conditions will have to be prescribed.


  #-----------------------
  #- nondimensionalization

  local adim_Sche = 1e0 / (P1 / R0 * dx * dy)
  ru_Sche .*= -adim_Sche
  rv_Sche .*= -adim_Sche

  ru2_Sche .= -999e0
  rv2_Sche .= -999e0

  for k in 1:NK
    for j in 1:NJ
      for i in 1:NI-1
        ru2_Sche[i, j, k] = ru_Sche[i+1, j, k]
      end
    end
  end
  for k in 1:NK
    for j in 1:NJ-1
      for i in 1:NI
        rv2_Sche[i, j, k] = rv_Sche[i, j+1, k]
      end
    end
  end

  #-----------------------
  #- boundary conditions

  ru2_Sche[0, :, :] = ru_Sche[1, :, :]
  ru2_Sche[NI, :, :] = ru_Sche[1, :, :]

  rv2_Sche[:, 0, :] .= 0e0
  rv2_Sche[:, NJ, :] .= 0e0

  #--------------------------------
  # ru2_Sche and rv2_Sche are now totally equivalent to grpifc and grpjfc.
  # Now, we have to compute drpx and drpy.
  # The relation between grpifc/grpjfc and drpx/drpy is exactly the same as the one in the Song scheme.



  for k in 1:NK
    for j in 1:NJ
      for i in 1:NI
        ru3_Sche[i, j, k] = 0.5e0 * (ru2_Sche[i, j, k] / gi[i, j, k, 1] + ru2_Sche[i-1, j, k] / gi[i-1, j, k, 1])
      end
    end
  end

  for k in 1:NK
    for j in 1:NJ
      for i in 1:NI
        rv3_Sche[i, j, k] = 0.5e0 * (rv2_Sche[i, j, k] / gj[i, j, k, 2] + rv2_Sche[i, j-1, k] / gj[i, j-1, k, 2])
      end
    end
  end

  for k in 1:NK
    for j in 1:NJ
      for i in 1:NI
        ru4_Sche[i, j, k] = (ru3_Sche[i, j, k] * ux[i, j] + rv3_Sche[i, j, k] * vx[i, j])
      end
    end
  end

  for k in 1:NK
    for j in 1:NJ
      for i in 1:NI
        rv4_Sche[i, j, k] = (ru3_Sche[i, j, k] * uy[i, j] + rv3_Sche[i, j, k] * vy[i, j])
      end
    end
  end



  #---------------------------------
  # ru2_Sche and rv2_Sche are the grpifc and grpjfc
  # ru4_Sche and rv4_Sche are the drpx and drpy
  #---------------------------------


  #---------------------------------
  # _ End of the code.
  #---------------------------------


  #---------------------------------
  # To illustrate how the code works on a simple case, 
  # Let us set:
  #  * rho=R0+ constant
  #  * h varies in y.
  #
  # The pressure terms can be decomposed into:
  #
  # (              Eq. P.2) only pressure terms from NK-1 to 1:                         0           0.737      0.737      0.737    ...
  # (part of Eq. P.1 + P.2) pressure terms from NK-1 to 1 + the rho term in NK face:    0.737       1.474      1.474      1.474    ...
  # (part of Eq. P.1 + P.2) pressure terms from NK-1 to 1 + the g.z_w term:           -55.75      -55.02     -55.02     -55.02     ...
  # (        Eq. P.1 + P.2) all pressure terms:				            -55.01      -54.28     -54.28     -54.28     ...
  #                                                                                                        
  # Then, a local additional term is added:                                                               
  #                                                                                                        
  # (part of Eq. R.1)       only this term:				 	      0.737       0.         0.         0.       ...
  #                                                                                                        
  # (Eq. P.1, P.2 and R.1) In the end, the addition of all these terms gives:         -54.276     -54.285    -54.285   -54.285     ...
  #---------------------------------



end