function gradup3df()

  #REAL(kind=rc_kind)  :: dufpdx, dvfpdy, dwfpdz
  #REAL(kind=rc_kind)  :: dufdx , dvfdy , dwfdz
  #REAL(kind=rc_kind)  :: udpdx , vdpdy , wdpdz
  #
  #REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1)  :: ptc
  #REAL(kind=rc_kind), dimension(    0:NI,    NJ  ,   NK  )  :: ptfx
  #REAL(kind=rc_kind), dimension(      NI,  0:NJ,     NK  )  :: ptfy
  #REAL(kind=rc_kind), dimension(      NI,    NJ,   0:NK  )  :: ptfz
  local ptc = OffsetArrays.zeros(rc_kind, 0:NI+1, 0:NJ+1, 0:NK+1)
  diag_ptc(ptc)
  intpol_pt(ptc, ptfx, ptfy, ptfz)

  for k in 1:NK
    for j in 1:NJ
      for i in 1:NI

        dufpdx = uf[i, j, k] * ptfx[i, j, k] - uf[i-1, j, k] * ptfx[i-1, j, k]
        dvfpdy = vf[i, j, k] * ptfy[i, j, k] - vf[i, j-1, k] * ptfy[i, j-1, k]
        dwfpdz = wf[i, j, k] * ptfz[i, j, k] - wf[i, j, k-1] * ptfz[i, j, k-1]


        dufdx = (uf[i, j, k] - uf[i-1, j, k])
        dvfdy = (vf[i, j, k] - vf[i, j-1, k])
        dwfdz = (wf[i, j, k] - wf[i, j, k-1])


        udpdx = 0.5 * (uf[i, j, k] + uf[i-1, j, k]) * (ptfx[i, j, k] - ptfx[i-1, j, k])
        vdpdy = 0.5 * (vf[i, j, k] + vf[i, j-1, k]) * (ptfy[i, j, k] - ptfy[i, j-1, k])
        wdpdz = 0.5 * (wf[i, j, k] + wf[i, j, k-1]) * (ptfz[i, j, k] - ptfz[i, j, k-1])


      end
    end
  end

end
