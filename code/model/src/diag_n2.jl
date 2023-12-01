function diag_n2()
  local DLinv = 1e0 / DL
  for k in 1:NK
    for j in 1:NJ
      for i in 1:NI

        @views @inbounds rdy = 0.5e0 * ((rho[i+1, j+2, k+1] - sum(rho[:, j+2, k+1]) / NI) - (rho[i+1, j, k+1] - sum(rho[:, j, k+1]) / NI)) * vy[i+1, j+1] / LEN
        local rdz = 0e0
        local rpdz = 0e0
        if (k == NK)
          rdz = (rho[i+1, j+1, k+1] - rho[i+1, j+1, k]) * wz[i+1, j+1, k+1] * DLinv
          @views @inbounds rpdz = (rho[i+1, j+1, k+1] - sum(rho[:, j+1, k+1]) / NI - (rho[i+1, j+1, k] - sum(rho[:, j+1, k]) / NI)) * wz[i+1, j+1, k+1] * DLinv
        elseif (k == 1)
          rdz = (rho[i+1, j+1, k+2] - rho[i+1, j+1, k+1]) * wz[i+1, j+1, k+1] * DLinv
          @views @inbounds rpdz = (rho[i+1, j+1, k+2] - sum(rho[:, j+1, k+2]) / NI - (rho[i+1, j+1, k+1] - sum(rho[:, j+1, k+1]) / NI)) * wz[i+1, j+1, k+1] * DLinv
        else
          rdz = 0.5e0 * (rho[i+1, j+1, k+2] - rho[i+1, j+1, k]) * wz[i+1, j+1, k+1] * DLinv
          @views @inbounds rpdz = 0.5e0 * (rho[i+1, j+1, k+2] - sum(rho[:, j+1, k+2]) / NI - (rho[i+1, j+1, k] - sum(rho[:, j+1, k]) / NI)) * wz[i+1, j+1, k+1] * DLinv
        end

        freqby[i+1, j+1, k+1] = (-gpr * 10e0 / R0) * rdy
        freqbz[i+1, j+1, k+1] = (-gpr * 10e0 / R0) * rpdz
        freqN2[i+1, j+1, k+1] = (-gpr * 10e0 / R0) * rdz

      end
    end
  end
end