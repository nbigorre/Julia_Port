function linerelax_segment(chunk, nxm, nym, nzm, cp, p, fn)
  local info = Ref(0)
  local rh = zeros(rc_kind, NK)
  local subd = zeros(rc_kind, NK)
  local supd = zeros(rc_kind, NK)
  local dia = zeros(rc_kind, NK)
  for j in chunk
    local istart = div(NI, 2)
    for ii in 1:nxm
      local i = ii + istart
      if (i > nxm)
        i -= nxm
      end
      local im1 = i - 1
      local ip1 = i + 1
      if (i == 1)
        im1 = nxm
      end
      if (i == nxm)
        ip1 = 1
      end
      for k in 1:nzm
        rh[k] = -(cp[2, i, j, k] * p[1+ip1, j+1, k+1]
                  + cp[3, i, j, k] * p[1+i, j+2, k+1]
                  + cp[4, i, j, k] * p[1+im1, j+1, k+1]
                  + cp[5, i, j, k] * p[1+i, j, k+1]
                  + cp[8, i, j, k] * p[1+im1, j+2, k+1]
                  + cp[9, i, j, k] * p[1+im1, j, k+1]
                  + cp[10, i, j, k] * p[1+ip1, j, k+1]
                  + cp[11, i, j, k] * p[1+ip1, j+2, k+1]
                  + cp[12, i, j, k] * p[1+im1, j+1, k]
                  + cp[13, i, j, k] * p[1+ip1, j+1, k]
                  + cp[14, i, j, k] * p[1+ip1, j+1, k+2]
                  + cp[15, i, j, k] * p[1+im1, j+1, k+2]
                  + cp[16, i, j, k] * p[1+i, j, k]
                  + cp[17, i, j, k] * p[1+i, j+2, k]
                  + cp[18, i, j, k] * p[1+i, j+2, k+2]
                  + cp[19, i, j, k] * p[1+i, j, k+2]
                  - fn[i, j, k])
        if (k != 1)
          subd[k] = cp[7, i, j, k]
        end
        if (k != nzm)
          supd[k] = cp[6, i, j, k]
        end
        dia[k] = cp[1, i, j, k]
      end
      rh[1] -= cp[7, i, j, 1] * p[i+1, j+1, 1]
      rh[nzm] -= cp[6, i, j, nzm] * p[i+1, j+1, nzm+2]

      info[] = dgtsl(nzm, subd, dia, supd, rh)
      #@ccall "./PSOM_LIB.so".dgtsl_(Ref(nzm)::Ref{Int}, pointer(subd)::Ptr{rc_kind}, pointer(dia)::Ptr{rc_kind}, pointer(supd)::Ptr{rc_kind}, pointer(rh)::Ptr{rc_kind}, info::Ref{Int})::Cvoid
      if (info[] != 0)
        println("error in linerelax")
        exit(1)
      end
      @inbounds p[i+1, j+1, 2:nzm+1] .= rh[1:nzm]
    end
  end
end


function linerelax(nxm, nym, nzm, cp, p, fn)
  
  #local chunks = Iterators.partition(1:nym, div(nym, Threads.nthreads()))
  #local tasks = map(chunks) do chunk
  #  Threads.@spawn linerelax_segment(chunk, nxm, nym, nzm, cp, p, fn)
  #end
  #wait.(tasks)
  linerelax_segment(1:nym, nxm, nym, nzm, cp, p, fn)

end