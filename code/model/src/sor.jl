function sor(nxm, nym, nzm, cp, p, fn)
      local rlx = 1.8e0
      for iter in 1:8
            for k in 1:nzm
                  for j in 1:nym
                        for i in nxm
                              pstar= ( cp[2,i,j,k]*p[i,j-1,k-1] 
                                            +cp[3,i,j,k ]*p[i,j+1,k]
                                            +cp[4,i,j,k ]*p[i-1,j,k]
                                            +cp[5,i,j,k ]*p[i,j-1,k]
                                            +cp[6,i,j,k ]*p[i,j,k+1]
                                            +cp[7,i,j,k ]*p[i,j,k-1]
                                            +cp[8,i,j,k ]*p[i-1,j+1,k] 
                                            +cp[9,i,j,k ]*p[i-1,j-1,k] 
                                            +cp[10,i,j,k]*p[i+1,j-1,k]
                                            +cp[11,i,j,k]*p[i+1,j+1,k]
                                            +cp[12,i,j,k]*p[i-1,j,k-1]
                                            +cp[13,i,j,k]*p[i+1,j,k-1]
                                            +cp[14,i,j,k]*p[i+1,j,k+1]
                                            +cp[15,i,j,k]*p[i-1,j,k+1]
                                            +cp[16,i,j,k]*p[i,j-1,k-1]
                                            +cp[17,i,j,k]*p[i,j+1,k-1]
                                            +cp[18,i,j,k]*p[i,j+1,k+1]
                                            +cp[19,i,j,k]*p[i,j-1,k+1]
                                            - fn[i,j,k] )/(-cp[1,i,j,k])                        
                                        p[i,j,k] = (1e0 -rlx)*p[i,j,k] +rlx*pstar 
                        end
                  end
            end
      end
end