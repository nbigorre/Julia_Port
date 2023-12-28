function efill(nxm, nym, nzm, p)
   for k in 1:nzm
      for i in 1:nxm
         p[i, 0, k] = p[i, 1, k]
         p[i, nym+1, k] = p[i, nym, k]
      end
   end

   for k in 1:nzm
      for j in 1:nym
         p[0, j, k] = p[1, j, k]
         p[nxm+1, j, k] = p[nxm, j, k]
      end
   end
   for j in 1:nym
      for i in 1:nxm
         p[i, j, 0] = p[i, j, 1]
      end
   end

   for j in 0:nym+1
      for i in 0:nxm+1
         p[i, j, nzm+1] = -p[i, j, nzm]
      end
   end

   for k in 1:nzm
      p[0, 0, k] = p[1, 1, k]
      p[nxm+1, 0, k] = p[nxm, 1, k]
      p[0, nym+1, k] = p[1, nym, k]
      p[nxm+1, nym+1, k] = p[nxm, nym, k]
   end

   for i in 1:nxm
      p[i, 0, 0] = p[i, 1, 1]
      p[i, nym+1, 0] = p[i, nym, 1]
   end

   for j in 1:nym
      p[0, j, 0] = p[1, j, 1]
      p[nxm+1, j, 0] = p[nxm, j, 1]
   end
   p[0, 0, 0] = p[1, 1, 1]
   p[0, nym+1, 0] = p[1, nym, 1]
   p[nxm+1, 0, 0] = p[nxm, 1, 1]
   p[nxm+1, nym+1, 0] = p[nxm, nym, 1]
end