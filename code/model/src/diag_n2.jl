function diag_n2()
  local DLinv = 1e0 / DL
  for k in 1:NK
    for j in 1:NJ
      for i in 1:NI

        @views @inbounds rdy = 0.5e0 * ((rho[i, j+1, k] - sum(rho[:, j+1, k]) / NI) - (rho[i, j-1, k] - sum(rho[:, j-1, k]) / NI)) * vy[i, j] / LEN
        local rdz = 0e0
        local rpdz = 0e0
        if (k == NK)
          rdz = (rho[i, j, k] - rho[i, j, k-1]) * wz[i, j, k] * DLinv
          @views @inbounds rpdz = (rho[i, j, k] - sum(rho[:, j, k]) / NI - (rho[i, j, k-1] - sum(rho[:, j, k-1]) / NI)) * wz[i, j, k] * DLinv
        elseif (k == 1)
          rdz = (rho[i, j, k+1] - rho[i, j, k]) * wz[i, j, k] * DLinv
          @views @inbounds rpdz = (rho[i, j, k+1] - sum(rho[:, j, k+1]) / NI - (rho[i, j, k] - sum(rho[:, j, k]) / NI)) * wz[i, j, k] * DLinv
        else
          rdz = 0.5e0 * (rho[i, j, k+1] - rho[i, j, k-1]) * wz[i, j, k] * DLinv
          @views @inbounds rpdz = 0.5e0 * (rho[i, j, k+1] - sum(rho[:, j, k+1]) / NI - (rho[i, j, k-1] - sum(rho[:, j, k-1]) / NI)) * wz[i, j, k] * DLinv
        end

        freqby[i, j, k] = (-gpr * 10e0 / R0) * rdy
        freqbz[i, j, k] = (-gpr * 10e0 / R0) * rpdz
        freqN2[i, j, k] = (-gpr * 10e0 / R0) * rdz

      end
    end
  end
end