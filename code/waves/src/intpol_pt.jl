function intpol_pt(varc, varfx, varfy, varfz)
   #!----------------------------------------------------                   
   #!     Use 3rd Order accurate interpolation                              
   #!     interpolate cx,cy,cz onto the cell faces : cxf                    
   #!     doesn't use 2nd order QUICK scheme of interpolation now           
   #!     n is the previous time level                                      
   #  USE header
   #      integer i,j,k,i0,im1,ip1,ip2,j0,jm1,jp1,jp2 
   #      REAL(kind=rc_kind) :: cx1,cx2,cy1,cy2,Jack,xa,xb 
   #
   #      REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1)  :: varc
   #      REAL(kind=rc_kind), dimension(    0:NI,    NJ  ,   NK  )  :: varfx
   #      REAL(kind=rc_kind), dimension(      NI,  0:NJ,     NK  )  :: varfy
   #      REAL(kind=rc_kind), dimension(      NI,    NJ,   0:NK  )  :: varfz
   #!                                                                       
   #!     we need to impose boundary conditions like uu_x+vu_y+wu_z=0       
   #!     at i=0 and NI solid boundaries. If we use just interpolation      
   #!     to get these at the faces, we will not get the right value.       
   #!                                                                       
   if ((NI < 4) || (NJ < 4))
      println("prob in intpol,insuff i pts")
      exit(1)
   end

   local xa = 9.0 / 16.0
   local xb = -1.0 / 16.0
   #     periodic ew boundaries                                            
   for i in 0:NI
      i0 = i
      im1 = i - 1
      ip1 = i + 1
      ip2 = i + 2
      if (i == 0)
         i0 = NI
         im1 = NI - 1
         ip1 = i + 1
         ip2 = i + 2
      elseif (i == 1)
         i0 = i
         im1 = NI
         ip1 = i + 1
         ip2 = i + 2
      elseif (i == NI)
         i0 = i
         im1 = i - 1
         ip1 = 1
         ip2 = 2
      elseif (i == (NI - 1))
         i0 = i
         im1 = i - 1
         ip1 = i + 1
         ip2 = 1
      else
         i0 = i
         im1 = i - 1
         ip1 = i + 1
         ip2 = i + 2
      end

      for k in 1:NK
         for j in 1:NJ
            #     do for face to the east of the i-th cell                          
            varfx[i, j, k] = (xa * (varc[i0, j, k] + varc[ip1, j, k]) + xb * (varc[im1, j, k] + varc[ip2, j, k])) + 0.5 * (varc[i0, j, k] + varc[ip1, j, k])
         end
      end
   end

   @static if (!cppdefs.periodicew)

      for k in 1:NK
         for j in 1:NJ
            #cxf(0,j,k) = ufbcw(j,k) 
            #cxf(NI,j,k)= ufbce(j,k) 
         end
      end
   end

   for j in 1:NJ-1
      j0 = j
      jm1 = j - 1
      jp1 = j + 1
      jp2 = j + 2
      xa = 9.0 / 16.0
      xb = -1.0 / 16.0
      if (j == 1)
         j0 = j
         jm1 = NJ
         jp1 = j + 1
         jp2 = j + 2
         xa = 0.5
         xb = 0.0
      elseif (j == NJ - 1)
         j0 = j
         jm1 = j - 1
         jp1 = NJ
         jp2 = 1
         xa = 0.5
         xb = 0.0
      else
         j0 = j
         jm1 = j - 1
         jp1 = j + 1
         jp2 = j + 2
         xa = 9.0 / 16.0
         xb = -1.0 / 16.0
      end

      for k in 1:NK
         for i in 1:NI
            #     do for face to the north of the j-th cell                         
            varfy[i, j, k] = xa * (varc[i, j0, k] + varc[i, jp1, k]) + xb * (varc[i, jm1, k] + varc[i, jp2, k])
         end
      end
   end

   for k in 1:NK
      for i in 1:NI
         #cyf(i,NJ,k)= vfbcn(i,k) 
         #cyf(i,0,k)= vfbcs(i,k) 
      end
   end
   for k in 1:NK-1
      for j in 1:NJ
         for i in 1:NI
            varfz[i, j, k] = 0.5 * (varc[i, j, k] + varc[i, j, k+1])
         end
      end
   end
   #     czf would be zero at k=0,NK impermeable boundaries                
end
