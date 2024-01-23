function correctbc()
      # !     ------------------------------------------------------------      
      # USE header
      #                                                                       
      #      integer i,j,k,nexit 
      #      REAL(kind=rc_kind) :: flowin,flowout,const 

      #     Outflow :  correct and update                                     
      #     ---------                                                         
      local flowin = 0e0
      local flowout = 0e0
      local nexit = 0
      for k in 1:NK
            for i in 1:NI
                  #            if (vfbcs[i,k] >= 0e0)                              
                  flowin = flowin + vfbcs[i, k]
                  #            else                                                       
                  #               flowout= flowout - vfbcs[i,k]                           
                  #               nexit= nexit +1                                         
                  #            end                                                      
            end
      end
      for k in NK-1:NK
            for i in 1:NI
                  #            if (vfbcn[i,k] <= 0e0)                               
                  #               flowin= flowin -vfbcn[i,k]                              
                  #            else                                                       
                  flowout = flowout + vfbcn[i, k]
                  nexit = nexit + 1
                  #            end                                 
            end
      end
      #                                                                       
      #     MODIFY OUTFLOW VELOCITIES - assumes rectilinear grid in computing 
      local constant = (flowin - flowout) / nexit
      #     modify only vfbcn in top 2 layers                                 

      for k in NK-1:NK
            for i in 1:NI
                  #            if (vfbcn[i,k] > 0e0)                               
                  vfbcn[i, k] = vfbcn[i, k] + constant
                  vnorth[i, k] = vfbcs[i, k] / (vy[i, NJ] * Jjfc[i, NJ, k])
                  #            end           
            end
      end
end