
function hsolve(h,oldh,hdt,dtime)
  local ch = zeros(rc_kind, (9, NI, NJ))
  local rhs = zeros(rc_kind, (NI, NJ))


  local tol::rc_kind = 1e-12
  local dti::rc_kind = 1e0 / dtime
  local kount::Int = 0

  chfine(dtime, ch, rhs)
  #@ccall "./PSOM_LIB.so".chfine_(Ref(dtime)::Ref{rc_kind}, pointer(ch)::Ptr{rc_kind}, pointer(rhs)::Ptr{rc_kind})::Cvoid

  local rlx::rc_kind = 1.72e0

  for iter in 1:100
    for j in 1:NJ
      local i = 1
      local hstar = ( ch[2, i, j] * h[i+2, j+1]
                    + ch[3, i, j] * h[NI+1, j+1] 
                    + ch[4, i, j] * h[i+1, j+2] 
                    + ch[5, i, j] * h[i+1, j] 
                    + ch[6, i, j] * h[i+2, j+2] 
                    + ch[7, i, j] * h[NI+1, j+2] 
                    + ch[8, i, j] * h[i+2, j] 
                    + ch[9, i, j] * h[NI+1, j] 
                    - rhs[i, j]) / (-ch[1, i, j])
      h[i+1, j+1] = (1e0 - rlx) * h[i+1, j+1] + rlx * hstar

      for i in 2:NI-1
        local hstar = ( ch[2, i, j] * h[i+2, j+1]
                      + ch[3, i, j] * h[i, j+1] 
                      + ch[4, i, j] * h[i+1, j+2] 
                      + ch[5, i, j] * h[i+1, j] 
                      + ch[6, i, j] * h[i+2, j+2] 
                      + ch[7, i, j] * h[i, j+2] 
                      + ch[8, i, j] * h[i+2, j] 
                      + ch[9, i, j] * h[i, j] 
                      - rhs[i, j]) / (-ch[1, i, j])
        h[i+1, j+1] = (1e0 - rlx) * h[i+1, j+1] + rlx * hstar
      end
      local i = NI
      local hstar = ( ch[2, i, j] * h[2, j+1]
                    + ch[3, i, j] * h[i, j+1] 
                    + ch[4, i, j] * h[i+1, j+2] 
                    + ch[5, i, j] * h[i+1, j] 
                    + ch[6, i, j] * h[2, j+2] 
                    + ch[7, i, j] * h[i, j+2] 
                    + ch[8, i, j] * h[2, j] 
                    + ch[9, i, j] * h[i, j] 
                    - rhs[i, j]) / (-ch[1, i, j])
      h[i+1, j+1] = (1e0 - rlx) * h[i+1, j+1] + rlx * hstar
    end

    for l in 1:3
      @ccall "./PSOM_LIB.so".hfill_(Ref(dtime)::Ref{rc_kind}, pointer(h)::Ptr{rc_kind})::Cvoid
    end

    local maxres = 0e0

    for j in 1:NJ
      local i = 1
      local res = rhs[i,j] - ( 
        ch[1, i, j] * h[i+1, j+1]
      + ch[2, i, j] * h[i+2, j+1]
      + ch[3, i, j] * h[NI+1, j+1] 
      + ch[4, i, j] * h[i+1, j+2] 
      + ch[5, i, j] * h[i+1, j] 
      + ch[6, i, j] * h[i+2, j+2] 
      + ch[7, i, j] * h[NI+1, j+2] 
      + ch[8, i, j] * h[i+2, j] 
      + ch[9, i, j] * h[NI+1, j])
      if (abs(res) > maxres)
        maxres = abs(res)
      end
      for i in 2:NI-1
        local res = rhs[i,j] - ( 
          ch[1, i, j] * h[i+1, j+1]
        + ch[2, i, j] * h[i+2, j+1]
        + ch[3, i, j] * h[i, j+1] 
        + ch[4, i, j] * h[i+1, j+2] 
        + ch[5, i, j] * h[i+1, j] 
        + ch[6, i, j] * h[i+2, j+2] 
        + ch[7, i, j] * h[i, j+2] 
        + ch[8, i, j] * h[i+2, j] 
        + ch[9, i, j] * h[i, j])
        if (abs(res) > maxres)
          maxres = abs(res)
        end
      end
      local i = NI
      local res = rhs[i,j] - ( 
        ch[1, i, j] * h[i+1, j+1]
      + ch[2, i, j] * h[2, j+1]
      + ch[3, i, j] * h[i, j+1] 
      + ch[4, i, j] * h[i+1, j+2] 
      + ch[5, i, j] * h[i+1, j] 
      + ch[6, i, j] * h[2, j+2] 
      + ch[7, i, j] * h[i, j+2] 
      + ch[8, i, j] * h[2, j] 
      + ch[9, i, j] * h[i, j])
      if (abs(res) > maxres)
        maxres = abs(res)
      end
      if (maxres > 3000)
        println("hsolve.jl: STOP. res too large, i,j,maxres=",i,",",j,",",maxres)
        exit(1)
      end
    end


    if (maxres < tol)
      break
    end
  end

  for l in 1:1
    @ccall "./PSOM_LIB.so".mprove_(pointer(h)::Ptr{rc_kind}, pointer(ch)::Ptr{rc_kind}, pointer(rhs)::Ptr{rc_kind}, Ref(dtime)::Ref{rc_kind})::Cvoid
  end

end