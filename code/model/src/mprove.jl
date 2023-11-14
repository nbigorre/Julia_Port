function mprove(h,ch,oldrhs,dtime)
  local he = zeros(rc_kind, (NI+2, NJ+2))
  local rhs = zeros(rc_kind, (NI, NJ))
  local tol = 1e-10
  he .= 0e0

  for j in 1:NJ
    local i = 1
    rhs[i, j] = (ch[1,i,j] * h[i+1, j+1]
               + ch[2,i,j] * h[i+2, j+1]
               + ch[3,i,j] * h[NI+1, j+1]
               + ch[4,i,j] * h[i+1, j+2]
               + ch[5,i,j] * h[i+1, j]
               + ch[6,i,j] * h[i+2, j+2]
               + ch[7,i,j] * h[NI+1, j+2]
               + ch[8,i,j] * h[i+2, j]
               + ch[9,i,j] * h[NI+1, j]) - oldrhs[i, j]
    for i in 2:NI-1
      rhs[i, j] = (ch[1,i,j] * h[i+1, j+1]
                 + ch[2,i,j] * h[i+2, j+1]
                 + ch[3,i,j] * h[i, j+1]
                 + ch[4,i,j] * h[i+1, j+2]
                 + ch[5,i,j] * h[i+1, j]
                 + ch[6,i,j] * h[i+2, j+2]
                 + ch[7,i,j] * h[i, j+2]
                 + ch[8,i,j] * h[i+2, j]
                 + ch[9,i,j] * h[i, j]) - oldrhs[i, j]
    end
    local i = NI
    rhs[i, j] = (ch[1,i,j] * h[i+1, j+1]
               + ch[2,i,j] * h[2, j+1]
               + ch[3,i,j] * h[i, j+1]
               + ch[4,i,j] * h[i+1, j+2]
               + ch[5,i,j] * h[i+1, j]
               + ch[6,i,j] * h[2, j+2]
               + ch[7,i,j] * h[i, j+2]
               + ch[8,i,j] * h[2, j]
               + ch[9,i,j] * h[i, j]) - oldrhs[i, j]
  end

  local rlx = 1.72e0
  for iter in 1:1000
    for j in 1:NJ
      local i = 1
      local hstar = (ch[2,i,j] * he[i+2, j+1]
                   + ch[3,i,j] * he[NI, j+1]
                   + ch[4,i,j] * he[i+1, j+2]
                   + ch[5,i,j] * he[i+1, j]
                   + ch[6,i,j] * he[i+2, j+2]
                   + ch[7,i,j] * he[NI, j+2]
                   + ch[8,i,j] * he[i+2, j]
                   + ch[9,i,j] * he[NI, j] - rhs[i, j]) / (-ch[1,i,j])
      he[i+1,j+1] = (1e0 - rlx) * he[i+1, j+1] + rlx * hstar
      for i in 2:NI-1
        local hstar = (ch[2,i,j] * he[i+2, j+1]
                     + ch[3,i,j] * he[i, j+1]
                     + ch[4,i,j] * he[i+1, j+2]
                     + ch[5,i,j] * he[i+1, j]
                     + ch[6,i,j] * he[i+2, j+2]
                     + ch[7,i,j] * he[i, j+2]
                     + ch[8,i,j] * he[i+2, j]
                     + ch[9,i,j] * he[i, j] - rhs[i, j]) / (-ch[1,i,j])
        he[i+1,j+1] = (1e0 - rlx) * he[i+1, j+1] + rlx * hstar
      end
      local i = NI
      local hstar = (ch[2,i,j] * he[2, j+1]
                     + ch[3,i,j] * he[i, j+1]
                     + ch[4,i,j] * he[i+1, j+2]
                     + ch[5,i,j] * he[i+1, j]
                     + ch[6,i,j] * he[2, j+2]
                     + ch[7,i,j] * he[i, j+2]
                     + ch[8,i,j] * he[2, j]
                     + ch[9,i,j] * he[i, j] - rhs[i, j]) / (-ch[1,i,j])
        he[i+1,j+1] = (1e0 - rlx) * he[i+1, j+1] + rlx * hstar
    end
    for i in 1:NI
      he[i+1, 1] = he[i+1, 2]
      he[i+1, NJ+2] = he[i+1, NJ+1]
    end
    for j in 0:NJ+1
      he[1, j+1] = he[NI+1, j+1]
      he[NI+2, j+1] = he[2, j+1]
    end

    local maxres = 0e0
    for j in 1:NJ
      local i = 1
      local res = (ch[1,i,j] * he[i+1, j+1]
               + ch[2,i,j] * he[i+2, j+1]
               + ch[3,i,j] * he[NI+1, j+1]
               + ch[4,i,j] * he[i+1, j+2]
               + ch[5,i,j] * he[i+1, j]
               + ch[6,i,j] * he[i+2, j+2]
               + ch[7,i,j] * he[NI+1, j+2]
               + ch[8,i,j] * he[i+2, j]
               + ch[9,i,j] * he[NI+1, j])
      maxres = max(maxres, abs(res))
      for i in 2:NI-1
        local res = (ch[1,i,j] * he[i+1, j+1]
               + ch[2,i,j] * he[i+2, j+1]
               + ch[3,i,j] * he[i, j+1]
               + ch[4,i,j] * he[i+1, j+2]
               + ch[5,i,j] * he[i+1, j]
               + ch[6,i,j] * he[i+2, j+2]
               + ch[7,i,j] * he[i, j+2]
               + ch[8,i,j] * he[i+2, j]
               + ch[9,i,j] * he[i, j])
        maxres = max(maxres, abs(res))
      end
      local i = NI
      local res = (ch[1,i,j] * he[i+1, j+1]
               + ch[2,i,j] * he[2, j+1]
               + ch[3,i,j] * he[i, j+1]
               + ch[4,i,j] * he[i+1, j+2]
               + ch[5,i,j] * he[i+1, j]
               + ch[6,i,j] * he[2, j+2]
               + ch[7,i,j] * he[i, j+2]
               + ch[8,i,j] * he[2, j]
               + ch[9,i,j] * he[i, j])
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