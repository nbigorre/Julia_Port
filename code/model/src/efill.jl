function efill(nxm, nym, nzm, p)
   for k in 1:nzm
      for i in 1:nxm
         p[i+1, 1, k+1] = p[i+1, 2, k+1]
         p[i+1, nym+2, k+1] = p[i+1, nym+1, k+1]
      end
   end

   for k in 1:nzm
      for j in 1:nym
         p[1, j+1, k+1] = p[2, j+1, k+1]
         p[nxm+2, j+1, k+1] = p[nxm+1, j+1, k+1]
      end
   end
   for j in 1:nym
      for i in 1:nxm
         p[i+1, j+1, 1] = p[i+1, j+1, 2]
      end
   end

   for j in 0:nym+1
      for i in 0:nxm+1
         p[i+1, j+1, nzm+2] = -p[i+1, j+1, nzm+1]
      end
   end

   for k in 1:nzm
      p[1, 1, k+1] = p[2, 2, k+1]
      p[nxm+2, 1, k+1] = p[nxm+1, 2, k+1]
      p[1, nym+2, k+1] = p[2, nym+1, k+1]
      p[nxm+2, nym+2, k+1] = p[nxm+1, nym+1, k+1]
   end

   for i in 1:nxm
      p[i+1, 1, 1] = p[i+1, 2, 2]
      p[i+1, nym+2, 1] = p[i+1, nym+1, 2]
   end

   for j in 1:nym
      p[1, j+1, 1] = p[2, j+1, 2]
      p[nxm+2, j+1, 1] = p[nxm+1, j+1, 2]
   end
   p[1, 1, 1] = p[2, 2, 2]
   p[1, nym+2, 1] = p[2, nym+1, 2]
   p[nxm+2, 1, 1] = p[nxm+1, 2, 2]
   p[nxm+2, nym+2, 1] = p[nxm+1, nym+1, 2]
end