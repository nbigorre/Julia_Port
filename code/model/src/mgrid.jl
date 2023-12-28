function mgrid(p, dtime, edt, cfcdiv)
     local nx = zeros(Int, ngrid)
     local ny = zeros(Int, ngrid)
     local nz = zeros(Int, ngrid)
     local ntint = zeros(Int, ngrid)
     local ntout = zeros(Int, ngrid)
     local loco = zeros(Int, ngrid)
     local loci = zeros(Int, ngrid)
     local loccp = zeros(Int, ngrid)

     local cp = zeros(rc_kind, 19 * maxint)
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
     ntout[1] = (nx[1] + 2) * (ny[1] + 2) * (nz[1] + 2)

     loco[1] = 1
     loci[1] = 1
     loccp[1] = 1
     for m in 2:ngrid
          nx[m] = div(nx[m-1], 2)
          ny[m] = div(ny[m-1], 2)
          nz[m] = div(nz[m-1], 2)
          ntint[m] = nx[m] * ny[m] * nz[m]
          ntout[m] = (nx[m] + 2) * (ny[m] + 2) * (nz[m] + 2)
          loco[m] = loco[m-1] + ntout[m-1]
          loci[m] = loci[m-1] + ntint[m-1]
          loccp[m] = loccp[m-1] + 19 * ntint[m-1]
     end
     

     cpfine(dtime, cp, rhs)
     #@ccall "./PSOM_LIB.so".cpfine_(Ref(dtime)::Ref{rc_kind}, pointer(cp)::Ptr{rc_kind}, pointer(rhs)::Ptr{rc_kind})::Cvoid

     
     for m in 2:ngrid
          local cpf = reshape(view(cp, loccp[m-1]:(loccp[m-1]-1+(19*nx[m]*2*ny[m]*2*nz[m]*2))), 19, nx[m] * 2, ny[m] * 2, nz[m] * 2)
          local cpc = reshape(view(cp, loccp[m]:(loccp[m]-1+(19*nx[m]*ny[m]*nz[m]))), 19, nx[m], ny[m], nz[m])
          cpcors(nx[m], ny[m], nz[m], cpf, cpc)
          #@ccall "./PSOM_LIB.so".cpcors_(Ref(nx[m])::Ref{Int}, Ref(ny[m])::Ref{Int}, Ref(nz[m])::Ref{Int}, pointer(cp,loccp[m-1])::Ptr{rc_kind}, pointer(cp,loccp[m])::Ptr{rc_kind})::Cvoid
     end
     

     for ncycle in 1:noc
          for l in loco[2]:maxout
               p[l] = 0e0
          end
          local m = 1

          for l in 1:nu1
               
               local cp_r = reshape(view(cp, loccp[m]:(loccp[m]-1+(19*nx[m]*ny[m]*nz[m]))), 19, nx[m], ny[m], nz[m])
               local p_r = OffsetArrays.reshape(view(p, loco[m]:(loco[m]-1+((nx[m]+2)*(ny[m]+2)*(nz[m]+2)))), 0:nx[m] + 1, 0:ny[m] + 1, 0:nz[m] + 1)
               local rhs_r = reshape(view(rhs, loci[m]:(loci[m]-1+(nx[m]*ny[m]*nz[m]))), nx[m], ny[m], nz[m])
               linerelax(nx[m], ny[m], nz[m], cp_r, p_r, rhs_r)
               #@ccall "./PSOM_LIB.so".linerelax_(Ref(nx[m])::Ref{Int}, Ref(ny[m])::Ref{Int}, Ref(nz[m])::Ref{Int}, pointer(cp,loccp[m])::Ptr{rc_kind}, pointer(p, loco[m])::Ptr{rc_kind}, pointer(rhs, loci[m])::Ptr{rc_kind})::Cvoid

               mgpfill(dtime, reshape(view(p, 1:((NI+2)*(NJ+2)*(NK+2))), (NI + 2), (NJ + 2), (NK + 2)))
               #@ccall "./PSOM_LIB.so".mgpfill_(Ref(dtime)::Ref{rc_kind}, pointer(p)::Ptr{rc_kind})::Cvoid
          end

          let
               local cp_r = reshape(view(cp, loccp[m]:(loccp[m]-1+(19*nx[m]*ny[m]*nz[m]))), 19, nx[m], ny[m], nz[m])
               local p_r = OffsetArrays.reshape(view(p, loco[m]:(loco[m]-1+((nx[m]+2)*(ny[m]+2)*(nz[m]+2)))), 0:nx[m] + 1, 0:ny[m] + 1, 0:nz[m] + 1)
               local rhs_r = reshape(view(rhs, loci[m]:(loci[m]-1+(nx[m]*ny[m]*nz[m]))), nx[m], ny[m], nz[m])
               local res_r = reshape(view(res, 1:(nx[m]*ny[m]*nz[m])), nx[m], ny[m], nz[m])
               maxres[] = resid(m, nx[m], ny[m], nz[m], cp_r, p_r, rhs_r, res_r)
               #@ccall "./PSOM_LIB.so".resid_(Ref(m)::Ref{Int}, Ref(nx[m])::Ref{Int}, Ref(ny[m])::Ref{Int}, Ref(nz[m])::Ref{Int}, pointer(cp, loccp[m])::Ptr{rc_kind}, pointer(p, loco[m])::Ptr{rc_kind}, pointer(rhs, loci[m])::Ptr{rc_kind}, pointer(res)::Ptr{rc_kind}, maxres::Ref{rc_kind})::Cvoid
          end

          if (maxres[] < tol)
               return
          end
          let
               local res_r = reshape(view(res, 1:(nx[m]*ny[m]*nz[m])), nx[m], ny[m], nz[m])
               local rhs_rp1 = reshape(view(rhs, loci[m+1]:(loci[m+1]-1+(div(nx[m],2)*div(ny[m],2)*div(nz[m],2)))), div(nx[m],2), div(ny[m],2), div(nz[m],2))
               restrict(nx[m], ny[m], nz[m], res_r, rhs_rp1)
               #@ccall "./PSOM_LIB.so".restrict_(Ref(nx[m])::Ref{Int}, Ref(ny[m])::Ref{Int}, Ref(nz[m])::Ref{Int}, pointer(res)::Ptr{rc_kind}, pointer(rhs, loci[m+1])::Ptr{rc_kind})::Cvoid
          end
          for m in 2:ngrid-1
               local cp_r = reshape(view(cp, loccp[m]:(loccp[m]-1+(19*nx[m]*ny[m]*nz[m]))), 19, nx[m], ny[m], nz[m])
               local p_r = OffsetArrays.reshape(view(p, loco[m]:(loco[m]-1+((nx[m]+2)*(ny[m]+2)*(nz[m]+2)))), 0:nx[m] + 1, 0:ny[m] + 1, 0:nz[m] + 1)
               local rhs_r = reshape(view(rhs, loci[m]:(loci[m]-1+(nx[m]*ny[m]*nz[m]))), nx[m], ny[m], nz[m])
               local res_r = reshape(view(res, 1:(nx[m]*ny[m]*nz[m])), nx[m], ny[m], nz[m])
               local rhs_rp1 = reshape(view(rhs, loci[m+1]:(loci[m+1]-1+(div(nx[m],2)*div(ny[m],2)*div(nz[m],2)))), div(nx[m],2), div(ny[m],2), div(nz[m],2))
               for _ in 1:nu1
                    linerelax(nx[m], ny[m], nz[m], cp_r, p_r, rhs_r)
                    #@ccall "./PSOM_LIB.so".linerelax_(Ref(nx[m])::Ref{Int}, Ref(ny[m])::Ref{Int}, Ref(nz[m])::Ref{Int}, pointer(cp, loccp[m])::Ptr{rc_kind}, pointer(p, loco[m])::Ptr{rc_kind}, pointer(rhs, loci[m])::Ptr{rc_kind})::Cvoid
               end
               maxres[] = resid(m, nx[m], ny[m], nz[m], cp_r, p_r, rhs_r, res_r)
               #@ccall "./PSOM_LIB.so".resid_(Ref(m)::Ref{Int}, Ref(nx[m])::Ref{Int}, Ref(ny[m])::Ref{Int}, Ref(nz[m])::Ref{Int}, pointer(cp, loccp[m])::Ptr{rc_kind}, pointer(p, loco[m])::Ptr{rc_kind}, pointer(rhs, loci[m])::Ptr{rc_kind}, pointer(res)::Ptr{rc_kind}, maxres::Ref{rc_kind})::Cvoid
               restrict(nx[m], ny[m], nz[m], res_r, rhs_rp1)
               #@ccall "./PSOM_LIB.so".restrict_(Ref(nx[m])::Ref{Int}, Ref(ny[m])::Ref{Int}, Ref(nz[m])::Ref{Int}, pointer(res)::Ptr{rc_kind}, pointer(rhs, loci[m+1])::Ptr{rc_kind})::Cvoid
          end


          for m in ngrid:-1:2
               local cp_r = reshape(view(cp, loccp[m]:(loccp[m]-1+(19*nx[m]*ny[m]*nz[m]))), 19, nx[m], ny[m], nz[m])
               local p_r = OffsetArrays.reshape(view(p, loco[m]:(loco[m]-1+((nx[m]+2)*(ny[m]+2)*(nz[m]+2)))), 0:nx[m] + 1, 0:ny[m] + 1, 0:nz[m] + 1)
               local rhs_r = reshape(view(rhs, loci[m]:(loci[m]-1+(nx[m]*ny[m]*nz[m]))), nx[m], ny[m], nz[m])
               if m == ngrid
                    sor(nx[m], ny[m], nz[m], cp_r, p_r, rhs_r)
                    #@ccall "./PSOM_LIB.so".sor_(Ref(nx[m])::Ref{Int}, Ref(ny[m])::Ref{Int}, Ref(nz[m])::Ref{Int}, pointer(cp, loccp[m])::Ptr{rc_kind}, pointer(p, loco[m])::Ptr{rc_kind}, pointer(rhs, loci[m])::Ptr{rc_kind})::Cvoid
               else
                    for _ in 1:nu2
                         linerelax(nx[m], ny[m], nz[m], cp_r, p_r, rhs_r)
                         #@ccall "./PSOM_LIB.so".linerelax_(Ref(nx[m])::Ref{Int}, Ref(ny[m])::Ref{Int}, Ref(nz[m])::Ref{Int}, pointer(cp, loccp[m])::Ptr{rc_kind}, pointer(p, loco[m])::Ptr{rc_kind}, pointer(rhs, loci[m])::Ptr{rc_kind})::Cvoid
                    end
               end
               efill(nx[m], ny[m], nz[m], p_r)
               #@ccall "./PSOM_LIB.so".efill_(Ref(nx[m])::Ref{Int}, Ref(ny[m])::Ref{Int}, Ref(nz[m])::Ref{Int}, pointer(p, loco[m])::Ptr{rc_kind})::Cvoid

               let
                    local p_fin = reshape(view(p, loco[m-1]:(loco[m-1]-1+((2*nx[m]+2)*(2*ny[m]+2)*(2*nz[m]+2)))), (2 * nx[m] + 2), (2 * ny[m] + 2), (2 * nz[m] + 2))
                    prolong(nx[m], ny[m], nz[m], p_r, p_fin)
                    #@ccall "./PSOM_LIB.so".prolong_(Ref(nx[m])::Ref{Int}, Ref(ny[m])::Ref{Int}, Ref(nz[m])::Ref{Int}, pointer(p, loco[m])::Ptr{rc_kind}, pointer(p, loco[m-1])::Ptr{rc_kind})::Cvoid
               end
          end
     end
     

     local m = 1
     let
          local cp_r = reshape(view(cp, loccp[m]:(loccp[m]-1+(19*nx[m]*ny[m]*nz[m]))), 19, nx[m], ny[m], nz[m])
          local p_r = OffsetArrays.reshape(view(p, loco[m]:(loco[m]-1+((nx[m]+2)*(ny[m]+2)*(nz[m]+2)))), 0:nx[m] + 1, 0:ny[m] + 1, 0:nz[m] + 1)
          local p_base = reshape(view(p, 1:((NI+2)*(NJ+2)*(NK+2))), (NI + 2), (NJ + 2), (NK + 2))
          for l in 1:nu1
               linerelax(nx[m], ny[m], nz[m], cp_r, p_r, rhs_r)
               #@ccall "./PSOM_LIB.so".linerelax_(Ref(nx[m])::Ref{Int}, Ref(ny[m])::Ref{Int}, Ref(nz[m])::Ref{Int}, pointer(cp, loccp[m])::Ptr{rc_kind}, pointer(p, loco[m])::Ptr{rc_kind}, pointer(rhs, loci[m])::Ptr{rc_kind})::Cvoid
               
               mgpfill(dtime, p_base)
               #@ccall "./PSOM_LIB.so".mgpfill_(Ref(dtime)::Ref{rc_kind}, pointer(p)::Ptr{rc_kind})::Cvoid
          end
          local rhs_r = reshape(view(rhs, loci[m]:(loci[m]-1+(nx[m]*ny[m]*nz[m]))), nx[m], ny[m], nz[m])
          local res_r = reshape(view(res, 1:(nx[m]*ny[m]*nz[m])), nx[m], ny[m], nz[m])
          maxres[] = resid(m, nx[m], ny[m], nz[m], cp_r, p_r, rhs_r, res_r)
          #@ccall "./PSOM_LIB.so".resid_(Ref(m)::Ref{Int}, Ref(nx[m])::Ref{Int}, Ref(ny[m])::Ref{Int}, Ref(nz[m])::Ref{Int}, pointer(cp, loccp[m])::Ptr{rc_kind}, pointer(p, loco[m])::Ptr{rc_kind}, pointer(rhs, loci[m])::Ptr{rc_kind}, pointer(res)::Ptr{rc_kind}, maxres::Ref{rc_kind})::Cvoid
     end
     maxres[] = maxres[] / edt

     
end