function ini_st()
  #
  #  USE header
  #
  #!  initializes TS profiles from the Gulf Stream
  #!  regrided for the model
  #  implicit none
  #  integer n,i,j,k
  #
  #! maximum and minimum potential density values
  #  real(kind=rc_kind), parameter :: rhoS=1025.8d0
  #
  #! front parameters
  #!    + Af :length of the front in the y and z directions
  #!    + PI: inflection points position of the tanh functions:
  #!      + yPI:horizontal position of the center of the front.
  #!      + zPI:same as yPI. Allows to vary the mld along y
  #!    + mld:mixed layer depth
  #
  #  real(kind=rc_kind) :: Afy, mld, Afz 
  #  real(kind=rc_kind) :: Arho, yPI, zPI, z, Fz, Fy(NJ)
  local rhoS=1025.8e0
  local yPI = 0.5 * (yc[NJ+1] + yc[0])
  local Fy = zeros(rc_kind, NJ)
  # SUBMESOSCALE front

  local Afy = 15e0 / 4e0
  local Afy = 40e0 / 4e0
  local mld = 50e0
  local mld = 50e0
  local Afz = 3e2 / 4e0

  local zPI = mld + Afz

  for j in 1:NJ
    #Fy(j)=(tanh( (yc(j)-yPI)/Afy  ))*1.25
    Fy[j] = (tanh((yc[j] - yPI) / Afy)) / (4e0 / 3e0)
    #print*, j, Fy(j)
  end

  for k in 0:NK+1
    for j in 1:NJ
      z = DL * zc[1, j, k]
      #Fz=(tanh( (zPI*Fy(j)-z)/Afz )+1e0)/4.
      Fz = (tanh((zPI * Fy[j] - z) / Afz) + 1e0) / 2e0
      #s(0,j,k,0)=s(0,j,k,0)+Fz
      s[0, j, k, 0] = rhoS + Fz
      #if (j==1 .or. j==NJ) print*, k, z, Fz
    end
  end

  for k in 1:NK
    s[0, 0, k, 0] = s[0, 1, k, 0]
    s[0, NJ+1, k, 0] = s[0, NJ, k, 0]
    #    do j=350,NJ+1
    #      s(0,j,k,0)=s(0,349,k,0)
    #    end
  end

  for k in 0:NK+1
    for j in 0:NJ+1
      for i in 0:NI+1
        s[i, j, k, 0] = s[0, j, k, 0]
      end
    end
  end

  # reference density values for sponge layers
  for k in 0:NK+1
    rho_refS[k] = s[0, 1, k, 0]
    rho_refN[k] = s[0, NJ, k, 0]
  end
  for j in 0:NJ+1
    rho_refB[j] = s[0, j, 0, 0]
  end

  return

end



