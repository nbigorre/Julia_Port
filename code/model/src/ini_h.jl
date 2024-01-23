function ini_h()
  #  !     ------------------------------------------------                  
  #  ! Initialize h using thermal wind balance
  #  ! Any density distribution is fine. 

  #  ! In case the bottom if flat, h is defined in such a way that it exactly cancels speed at the bottom, 
  #  !                             supposing thermal wind balace.
  #  ! otherwise, velocity will be canceled at a given geopotential depth.

  local sig = Ref(1)
  local sigup = Ref(1e0)
  local fu = OffsetArrays.zeros(rc_kind, 0:NJ)
  if (lv_flat_bottom == 1)
    # Two strategies are possible.
    # optionl=1 works for solid boundaries in y=0 and y=ymax
    # optionl=2 is based on the computation of the baroclinic pressure term. It cancels velocity at the CENTER
    #                    of the lowermost cell.
    local optionl = 1

    if (optionl == 1)
      local constant = DL / R0 #fu is actually fu / g and is at cell faces
      for i in 1:NI
        for j in 1:NJ-1
          fu[j] = 0e0
          for k in 1:NK
            # findzall has been called, so zf can be used. sigma not yet called, so wz not usable
            local dz = zf[i, j, k] - zf[i, j, k-1]
            local drho = (rho[i, j+1, k] - rho[i, j, k]) * vy[i, j] / LEN
            @fortSet("drho", drho, rc_kind)
            fu[j] = fu[j] + constant * drho * dz  #fu is at cell faces,is actually fu/g
          end
        end
        fu[0] = fu[1]
        fu[NJ] = fu[NJ-1]
        # at k=NK, fu = g*hy
        h[i, 0] = 0e0
        for j in 1:NJ+1
          h[i, j] = h[i, j-1] - (LEN / vy[i, j]) * fu[j-1]
        end
      end
    end

    if (optionl == 2)
      for i in 1:NI+1
        h[i, 0] = 0e0
        for j in 1:NJ
          h[i, j] = h[i, j-1] - HL / (gpr * gj[i, j-1, 1, 2]) * grpjfc[i, j-1, 1]
        end
      end
    end

  else

    # In the general case of non-flat bottom, h cannot cancel bottom speed.
    # Here, it simply cancels speed at a given depth, zl, that has to be within the domain at every (x,y) location.
    # It only works for solid boundaries in y=0 and y=ymax and in the case where h is flat at y=0.

    staticsigma()
    #@ccall "./PSOM_LIB.so".staticsigma_()::Cvoid
    rpevalgrad_Sche(0)
    #@ccall "./PSOM_LIB.so".rpevalgrad_sche_(Ref(0)::Ref{Int})::Cvoid

    optionl = 1

    if (optionl == 1)
      local zl = -97e0
      for i in 1:NI
        h[i, 0] = 0e0
        h[i, 1] = 0e0
        for j in 2:NJ
          (sig[], sigup[]) = findsigma(i, j-1, zl)
          #@ccall "./PSOM_LIB.so".findsigma_(Ref(i)::Ref{Int}, Ref(j - 1)::Ref{Int}, Ref(zl)::Ref{rc_kind}, sig::Ref{Int}, sigup::Ref{rc_kind})::Cvoid# sig=1;sigup=1. 
          h[i, j] = h[i, j-2] - 1e0 / (0.5 * gpr * vy[i, j-1]) * (rv4_Sche[i, j-1, sig[]] * sigup[] + rv4_Sche[i, j-1, sig[]-1] * (1e0 - sigup[]))
        end
      end
    else
      println("optionl=1 is the only option in inith!")
    end

  end

  # ------------------------------------------------

  # Additionally, a barotropic jet can be added:
  #@ccall "./PSOM_LIB.so".ini_h_()::Cvoid


  # periodicity
  h[:, NJ+1] .= h[:, NJ]
  for j in 0:NJ+1
    h[0, j] = h[NI, j]
    h[NI+1, j] = h[1, j]
  end

  # h should be zero-average.
  local hsum = 0e0
  for i in 1:NI
    for j in 1:NJ
      hsum = hsum + h[i, j]
    end
  end
  global hmean = hsum / rc_kind(NI * NJ)
  @views @. h -= hmean
  @views @. h /= HL
  

  # h=0.

end
