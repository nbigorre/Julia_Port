function smooth()                                                             
      for j in 1:NJ
            for i in 1:NI
                  local Ddu = 0.5e0 * (D[i+1, j] - D[i-1, j])
                  local Ddv = 0.5e0 * (D[i, j+1] - D[i, j-1])
                  Ddx[i, j] = Ddu * ux[i, j] + Ddv * vx[i, j]
                  Ddy[i, j] = Ddu * uy[i, j] + Ddv * vy[i, j]
            end
      end
      local i = 0
      for j in 1:NJ
            local Ddu = D[i+1, j] - D[i, j]
            local Ddv = 0.5e0 * (D[i, j+1] - D[i, j-1])
            Ddx[i, j] = Ddu * ux[i, j] + Ddv * vx[i, j]
            Ddy[i, j] = Ddu * uy[i, j] + Ddv * vy[i, j]
      end
      i = NI + 1
      for j in 1:NJ
            local Ddu = D[i, j] - D[i-1, j]
            local Ddv = 0.5e0 * (D[i, j+1] - D[i, j-1])
            Ddx[i, j] = Ddu * ux[i, j] + Ddv * vx[i, j]
            Ddy[i, j] = Ddu * uy[i, j] + Ddv * vy[i, j]
      end
      local j = 0
      for i in 1:NI
            local Ddu = 0.5e0 * (D[i+1, j] - D[i-1, j])
            local Ddv = D[i, j+1] - D[i, j]
            Ddx[i, j] = Ddu * ux[i, j] + Ddv * vx[i, j]
            Ddy[i, j] = Ddu * uy[i, j] + Ddv * vy[i, j]
      end
      j = NJ + 1
      for i in 1:NI
            local Ddu = 0.5e0 * (D[i+1, j] - D[i-1, j])
            local Ddv = D[i, j] - D[i, j-1]
            Ddx[i, j] = Ddu * ux[i, j] + Ddv * vx[i, j]
            Ddy[i, j] = Ddu * uy[i, j] + Ddv * vy[i, j]
      end

      Ddx[0, 0] = 0.5e0 * (Ddx[1, 0] + Ddx[0, 1])
      Ddy[0, 0] = 0.5e0 * (Ddy[1, 0] + Ddy[0, 1])
      Ddx[NI+1, 0] = 0.5e0 * (Ddx[NI, 0] + Ddx[NI+1, 1])
      Ddy[NI+1, 0] = 0.5e0 * (Ddy[NI, 0] + Ddy[NI+1, 1])
      Ddx[NI+1, NJ+1] = 0.5e0 * (Ddx[NI, NJ+1] + Ddx[NI+1, NJ])
      Ddy[NI+1, NJ+1] = 0.5e0 * (Ddy[NI, NJ+1] + Ddy[NI+1, NJ])
      Ddx[0, NJ+1] = 0.5e0 * (Ddx[1, NJ+1] + Ddx[0, NJ])
      Ddy[0, NJ+1] = 0.5e0 * (Ddy[1, NJ+1] + Ddy[0, NJ])

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
