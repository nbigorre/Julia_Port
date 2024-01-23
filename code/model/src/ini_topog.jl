function ini_topog()
   #  USE header
   #  !     dy is the y grid spacing in m                                     
   #  !     creates the bottom topography                                     
   #  !     fits a spline to the described topog and computed Ddx(i)          
   #      integer i,j,NJnew,NJ2,jj,nfac,nfacm1,jtemp,inc 
   #      parameter (nfac=4) 
   #      parameter (nfacm1= nfac+1) 
   #      parameter (NJnew=nfac*NJ) 
   #      parameter (NJ2= NJnew+2*nfac) 
   #      REAL(kind=rc_kind) :: dynew,dep,seval,LbyD
   #      REAL(kind=rc_kind) :: ynew(-nfacm1:NJnew+nfac),                        &
   #     &     dnew(-nfacm1:NJnew+nfac),bd1(-nfacm1:NJnew+nfac),            &
   #     &     cd1(-nfacm1:NJnew+nfac),dd1(-nfacm1:NJnew+nfac)              
   #      REAL(kind=rc_kind) :: y,ytemp,deriv
   local nfac::Int = 4
   local nfacm1::Int = nfac+1
   local NJnew::Int = nfac*NJ
   local NJ2 = NJnew+2*nfac
   local dynew = dy / float(nfac)
   local LbyD = LEN / DL
   local DLinv = 1e0 / DL

   local ynew = OffsetArrays.zeros(rc_kind, -nfacm1:NJnew+nfac)
   local dnew = OffsetArrays.zeros(rc_kind, -nfacm1:NJnew+nfac)
   local bd1 = OffsetArrays.zeros(rc_kind, -nfacm1:NJnew+nfac)
   local cd1 = OffsetArrays.zeros(rc_kind, -nfacm1:NJnew+nfac)
   local dd1 = OffsetArrays.zeros(rc_kind, -nfacm1:NJnew+nfac)

   ynew[-nfacm1] = -(0.5 + float(nfacm1)) * dynew * 1e-3

   for j in -nfacm1+1:NJnew+nfac
      ynew[j] = ynew[j-1] + dynew * 1e-3
   end

   for j in -nfac:NJnew+nfac
      local dep = 0e0
      if (ynew[j] < 110)
         dep = 101e0    # Flat shelf 
      # dep= 101.+2*(ynew(j)-110.) ! Very close to Gawarkiewicz profile 
      elseif (ynew[j] < 118e0)
         dep = 101e0 + 1.5 * (ynew[j] - 108e0)^2 - 0.1 * (ynew[j] - 118e0)^2
      elseif (ynew[j] < 140e0)
         dep = 310e0 + 30 * (ynew[j] - 120e0)
      elseif (ynew[j] < 150e0)
         dep = 1060e0 - 1.5 * (ynew[j] - 150e0)^2
      else
         dep = 1060e0
      end

      dnew[j] = dep
   end

   #spline(NJ2, ynew, dnew, bd1, cd1, dd1)
   @ccall "./PSOM_LIB.so".spline_(Ref(NJ2)::Ref{Int}, pointer(ynew)::Ptr{rc_kind}, pointer(dnew)::Ptr{rc_kind}, pointer(bd1)::Ptr{rc_kind}, pointer(cd1)::Ptr{rc_kind}, pointer(dd1)::Ptr{rc_kind})::Cvoid

   for j in 0:NJ+1
      local dep = seval(NJ2, yc[j], ynew, dnew, bd1, cd1, dd1)
      #local dep = @ccall "./PSOM_LIB.so".seval_(Ref(NJ2)::Ref{Int}, Ref(yc[j])::Ref{rc_kind}, pointer(ynew)::Ptr{rc_kind}, pointer(dnew)::Ptr{rc_kind}, pointer(bd1)::Ptr{rc_kind}, pointer(cd1)::Ptr{rc_kind}, pointer(dd1)::Ptr{rc_kind})::rc_kind

      for jj in -nfacm1:NJnew+nfac
         if ((yc[j] < ynew[jj+1]) && (yc[j] >= ynew[jj]))
            local ytemp = yc[j] - ynew[jj]
            local jtemp::Int = nfac * j
            local inc::Int = div(nfac, 2)
            local deriv = 0e0
            if (j != NJ + 1)
               deriv = (dnew[jtemp+inc] - dnew[jtemp-inc]) / (ynew[jtemp+inc] - ynew[jtemp-inc])
            else
               deriv = (dnew[jtemp] - dnew[jtemp-inc]) / (ynew[jtemp] - ynew[jtemp-inc])
            end
            for i in 0:NI+1
               D[i, j] = -dep * DLinv
               Ddy[i, j] = deriv * LbyD * 1e-3
               Ddx[i, j] = 0e0
               #     multiplication by (L/D)*1.d-3 takes care of converting from km to 
               #     and non-dimensionalizing                                         
            end
            @goto loop_valid
         end
      end
      println("some prob in ini_topog")
      println(yc[j], jj)
      exit(1)
      @label loop_valid
   end

   @static if (cppdefs.file_output)
      local fs = open(string(dirout, "/Ddxdy.dat"), "w")
      for j in 0:NJ+1
         write(fs, string("j=    ", j, "\n"))
         for i in 0:NI+1
            write(fs, string(D[i, j], "  ", Ddx[i, j], "  ", Ddy[i, j], "\n"))
         end
      end
      close(fs)
   end

end