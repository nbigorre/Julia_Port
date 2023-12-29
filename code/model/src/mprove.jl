function mprove(h,ch,oldrhs,dtime)
  local he = OffsetArray(zeros(rc_kind, (NI+2, NJ+2)), 0:NI+1, 0:NJ+1)
  local rhs = zeros(rc_kind, (NI, NJ))
  local tol = 1e-10
  he .= 0e0

  for j in 1:NJ
    local i = 1
    rhs[i, j] = (ch[1,i,j] * h[i, j]
               + ch[2,i,j] * h[i+1, j]
               + ch[3,i,j] * h[NI, j]
               + ch[4,i,j] * h[i, j+1]
               + ch[5,i,j] * h[i, j-1]
               + ch[6,i,j] * h[i+1, j+1]
               + ch[7,i,j] * h[NI, j+1]
               + ch[8,i,j] * h[i+1, j-1]
               + ch[9,i,j] * h[NI, j-1]) - oldrhs[i, j]
    for i in 2:NI-1
      rhs[i, j] = (ch[1,i,j] * h[i, j]
                 + ch[2,i,j] * h[i+1, j]
                 + ch[3,i,j] * h[i-1, j]
                 + ch[4,i,j] * h[i, j+1]
                 + ch[5,i,j] * h[i, j-1]
                 + ch[6,i,j] * h[i+1, j+1]
                 + ch[7,i,j] * h[i-1, j+1]
                 + ch[8,i,j] * h[i+1, j-1]
                 + ch[9,i,j] * h[i-1, j-1]) - oldrhs[i, j]
    end
    local i = NI
    rhs[i, j] = (ch[1,i,j] * h[i, j]
               + ch[2,i,j] * h[1, j]
               + ch[3,i,j] * h[i-1, j]
               + ch[4,i,j] * h[i, j+1]
               + ch[5,i,j] * h[i, j-1]
               + ch[6,i,j] * h[1, j+1]
               + ch[7,i,j] * h[i-1, j+1]
               + ch[8,i,j] * h[1, j-1]
               + ch[9,i,j] * h[i-1, j-1]) - oldrhs[i, j]
  end

  local rlx = 1.72e0
  for iter in 1:1000
    for j in 1:NJ
      local i = 1
      local hstar = (ch[2,i,j] * he[i+1, j]
                   + ch[3,i,j] * he[NI, j]
                   + ch[4,i,j] * he[i, j+1]
                   + ch[5,i,j] * he[i, j-1]
                   + ch[6,i,j] * he[i+1, j+1]
                   + ch[7,i,j] * he[NI, j+1]
                   + ch[8,i,j] * he[i+1, j-1]
                   + ch[9,i,j] * he[NI, j-1] - rhs[i, j]) / (-ch[1,i,j])
      he[i,j] = (1e0 - rlx) * he[i, j] + rlx * hstar
      for i in 2:NI-1
        local hstar = (ch[2,i,j] * he[i+1, j]
                     + ch[3,i,j] * he[i-1, j]
                     + ch[4,i,j] * he[i, j+1]
                     + ch[5,i,j] * he[i, j-1]
                     + ch[6,i,j] * he[i+1, j+1]
                     + ch[7,i,j] * he[i-1, j+1]
                     + ch[8,i,j] * he[i+1, j-1]
                     + ch[9,i,j] * he[i-1, j-1] - rhs[i, j]) / (-ch[1,i,j])
        he[i,j] = (1e0 - rlx) * he[i, j] + rlx * hstar
      end
      local i = NI
      local hstar = (ch[2,i,j] * he[1, j]
                     + ch[3,i,j] * he[i-1, j]
                     + ch[4,i,j] * he[i, j+1]
                     + ch[5,i,j] * he[i, j-1]
                     + ch[6,i,j] * he[1, j+1]
                     + ch[7,i,j] * he[i-1, j+1]
                     + ch[8,i,j] * he[1, j-1]
                     + ch[9,i,j] * he[i-1, j-1] - rhs[i, j]) / (-ch[1,i,j])
        he[i,j] = (1e0 - rlx) * he[i, j] + rlx * hstar
    end
    for i in 1:NI
      he[i, 0] = he[i, 1]
      he[i, NJ+1] = he[i, NJ]
    end
    for j in 0:NJ+1
      he[0, j] = he[NI, j]
      he[NI+1, j] = he[1, j]
    end

    local maxres = 0e0
    for j in 1:NJ
      local i = 1
      local res = (ch[1,i,j] * he[i, j]
               + ch[2,i,j] * he[i+1, j]
               + ch[3,i,j] * he[NI, j]
               + ch[4,i,j] * he[i, j+1]
               + ch[5,i,j] * he[i, j-1]
               + ch[6,i,j] * he[i+1, j+1]
               + ch[7,i,j] * he[NI, j+1]
               + ch[8,i,j] * he[i+1, j-1]
               + ch[9,i,j] * he[NI, j-1])
      maxres = max(maxres, abs(res))
      for i in 2:NI-1
        local res = (ch[1,i,j] * he[i, j]
               + ch[2,i,j] * he[i+1, j]
               + ch[3,i,j] * he[i-1, j]
               + ch[4,i,j] * he[i, j+1]
               + ch[5,i,j] * he[i, j-1]
               + ch[6,i,j] * he[i+1, j+1]
               + ch[7,i,j] * he[i-1, j+1]
               + ch[8,i,j] * he[i+1, j-1]
               + ch[9,i,j] * he[i-1, j-1])
        maxres = max(maxres, abs(res))
      end
      local i = NI
      local res = (ch[1,i,j] * he[i, j]
               + ch[2,i,j] * he[1, j]
               + ch[3,i,j] * he[i-1, j]
               + ch[4,i,j] * he[i, j+1]
               + ch[5,i,j] * he[i, j-1]
               + ch[6,i,j] * he[1, j+1]
               + ch[7,i,j] * he[i-1, j+1]
               + ch[8,i,j] * he[1, j-1]
               + ch[9,i,j] * he[i-1, j-1])
        maxres = max(maxres, abs(res))
    end

    if (maxres > 3000)
      println("mprove.jl: STOP. res too large, i,j,maxres=",i,",",j,",",maxres)
      exit(1)
    end
    if (maxres < tol)
      break
    end
  end
  h .-= he
end