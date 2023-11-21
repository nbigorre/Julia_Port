function mgrid(p,dtime,edt,cfcdiv)
     local nx = zeros(Int, ngrid)
     local ny = zeros(Int, ngrid)
     local nz = zeros(Int, ngrid)
     local ntint = zeros(Int, ngrid)
     local ntout = zeros(Int, ngrid)
     local loco = zeros(Int, ngrid)
     local loci = zeros(Int, ngrid)
     local loccp = zeros(Int, ngrid)

     local cp = zeros(rc_kind, 19*maxint)
     local rhs = zeros(rc_kind, maxint)
     local res = zeros(rc_kind, int1)

     local noc = 30
     local nu1 = 2
     local nu2 = 1
     local tol_r = 1e-10

     nx[1] = NI
     ny[1] = NJ
     nz[1] = NK
     local tol = tol_r
     local maxres = Ref(0e0)
     local ncycle = 0
     if (cfcdiv <= tol)
          for m in 1:maxout
               p[m] = 0e0
          end
          maxres[] = edt * cfcdiv
          ncycle = 0
          return
     end

     tol = edt * tol
     ntint[1] = nx[1] * ny[1] * nz[1]
     ntout[1] = (nx[1]+2) * (ny[1]+2) * (nz[1]+2)

     loco[1] = 1
     loci[1] = 1
     loccp[1] = 1
     for m in 2:ngrid
          nx[m] = div(nx[m-1], 2)
          ny[m] = div(ny[m-1], 2)
          nz[m] = div(nz[m-1], 2)
          ntint[m] = nx[m] * ny[m] * nz[m]
          ntout[m] = (nx[m]+2) * (ny[m]+2) * (nz[m]+2)
          loco[m] = loco[m-1] + ntout[m - 1]
          loci[m] = loci[m-1] + ntint[m - 1]
          loccp[m] = loccp[m-1] + 19 * ntint[m - 1]
     end

     cpfine(dtime, cp, rhs)
     #@ccall "./PSOM_LIB.so".cpfine_(Ref(dtime)::Ref{rc_kind}, pointer(cp)::Ptr{rc_kind}, pointer(rhs)::Ptr{rc_kind})::Cvoid
     
     for m in 2:ngrid
          @ccall "./PSOM_LIB.so".cpcors_(Ref(nx[m])::Ref{Int}, Ref(ny[m])::Ref{Int}, Ref(nz[m])::Ref{Int}, pointer(cp,loccp[m-1])::Ptr{rc_kind}, pointer(cp,loccp[m])::Ptr{rc_kind})::Cvoid
     end

     for ncycle in 1:noc
          for l in loco[2]:maxout
               p[l] = 0e0
          end
          local m = 1

          for l in 1:nu1
               @ccall "./PSOM_LIB.so".linerelax_(Ref(nx[m])::Ref{Int}, pointer(ny,m)::Ptr{Int}, pointer(nz,m)::Ptr{Int}, pointer(cp,loccp[m])::Ptr{rc_kind}, pointer(p,loco[m])::Ptr{rc_kind}, pointer(rhs,loci[m])::Ptr{rc_kind})::Cvoid
               @ccall "./PSOM_LIB.so".mgpfill_(Ref(dtime)::Ref{rc_kind}, pointer(p)::Ptr{rc_kind})::Cvoid
          end
          @ccall "./PSOM_LIB.so".resid_(Ref(m)::Ref{Int}, Ref(nx[m])::Ref{Int}, Ref(ny[m])::Ref{Int}, Ref(nz[m])::Ref{Int}, pointer(cp,loccp[m])::Ptr{rc_kind}, pointer(p,loco[m])::Ptr{rc_kind}, pointer(rhs,loci[m])::Ptr{rc_kind}, pointer(res)::Ptr{rc_kind}, maxres::Ref{rc_kind})::Cvoid
          if (maxres[] < tol)
               return
          end
          @ccall "./PSOM_LIB.so".restrict_(Ref(nx[m])::Ref{Int}, Ref(ny[m])::Ref{Int}, Ref(nz[m])::Ref{Int}, pointer(res)::Ptr{rc_kind}, pointer(rhs,loci[m+1])::Ptr{rc_kind})::Cvoid

          for m in ngrid:-1:2
               if m == ngrid
                    @ccall "./PSOM_LIB.so".sor_(Ref(nx[m])::Ref{Int}, Ref(ny[m])::Ref{Int}, Ref(nz[m])::Ref{Int}, pointer(cp,loccp[m])::Ptr{rc_kind}, pointer(p,loco[m])::Ptr{rc_kind}, pointer(rhs,loci[m])::Ptr{rc_kind})::Cvoid
               else
                    for l in 1:nu2
                         @ccall "./PSOM_LIB.so".linerelax_(Ref(nx[m])::Ref{Int}, Ref(ny[m])::Ref{Int}, Ref(nz[m])::Ref{Int}, pointer(cp,loccp[m])::Ptr{rc_kind}, pointer(p,loco[m])::Ptr{rc_kind}, pointer(rhs,loci[m])::Ptr{rc_kind})::Cvoid
                    end
               end
               @ccall "./PSOM_LIB.so".efill_(Ref(nx[m])::Ref{Int}, Ref(ny[m])::Ref{Int}, Ref(nz[m])::Ref{Int}, pointer(p,loco[m])::Ptr{rc_kind})::Cvoid
               @ccall "./PSOM_LIB.so".prolong_(Ref(nx[m])::Ref{Int}, Ref(ny[m])::Ref{Int}, Ref(nz[m])::Ref{Int}, pointer(p,loco[m])::Ptr{rc_kind}, pointer(p,loco[m-1])::Ptr{rc_kind})::Cvoid
          end
     end

     local m = 1
     for l in 1:nu1
          @ccall "./PSOM_LIB.so".linerelax_(Ref(nx[m])::Ref{Int}, Ref(ny[m])::Ref{Int}, Ref(nz[m])::Ref{Int}, pointer(cp,loccp[m])::Ptr{rc_kind}, pointer(p,loco[m])::Ptr{rc_kind}, pointer(rhs,loci[m])::Ptr{rc_kind})::Cvoid
          @ccall "./PSOM_LIB.so".mgpfill_(Ref(dtime)::Ref{rc_kind}, pointer(p)::Ptr{rc_kind})::Cvoid
     end

     local m = 1
     @ccall "./PSOM_LIB.so".resid_(Ref(m)::Ref{Int}, Ref(nx[m])::Ref{Int}, Ref(ny[m])::Ref{Int}, Ref(nz[m])::Ref{Int}, pointer(cp,loccp[m])::Ptr{rc_kind}, pointer(p,loco[m])::Ptr{rc_kind}, pointer(rhs,loci[m])::Ptr{rc_kind}, pointer(res)::Ptr{rc_kind}, maxres::Ref{rc_kind})::Cvoid

     maxres[] = maxres[] / edt

end