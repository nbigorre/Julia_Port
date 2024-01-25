function sigma2z(x, y, sigma)
   #  !     -------------------------------------------------------           
   #  USE header, only: D,HDL,dztop,NK,h,pfac,rc_kind
   #  !     finds the value of z (non-dim by DL) given the value of sigma.    
   #  !     ht and dep are non-dim by DL.                                     
   #  !     At every (i,j), the column is divided into NK equal-depth cells.  
   #  REAL(kind=rc_kind) :: sigma,dnkm1,hpd,xfac,x,y,di,dj,h1,h2,z                     
   #  integer :: i,j
   #  !    pfac is the stretching in z. higher pfac gives more points near surf.
   #
   #  !      pfac= 5.d0  !c=NK32,dztop=0.5m                                   
   #  !==      pfac= 4.d0  !c=NK32,dztop=0.5m USED FOR SEVERAL MLI runs       
   #  !=      pfac=3.d0  !c used with NK=32 for first set of runs             
   local pfac = 2.0e0
   local dnkm1 = rc_kind(NK - 1)

   local i = Int(x)
   local j = Int(y)
   local di = x - i
   local dj = y - j
   if (sigma >= dnkm1)
      #In the surface layer
      local h1 = (h[i+1, j] - h[i, j]) * di + h[i, j]
      local h2 = (h[i+1, j+1] - h[i, j+1]) * di + h[i, j+1]
      local hpd = (h2 - h1) * dj + h1
      local hpd = hpd * HDL + dztop
      return (sigma - dnkm1) * hpd - dztop
   else
      #below the surface layer
      xfac = (dnkm1 - sigma) / dnkm1
      return (exp(pfac * xfac) - 1e0) * (D[i, j] + dztop) / (exp(pfac) - 1d0) - dztop
   end

   return 909090e0
end
