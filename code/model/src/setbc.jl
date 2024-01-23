function setbc(step)
      #     ------------------------------------------------                  
      # USE header
      #     sets the initial face velocities  and fills the boundary arrays.  
      #                                                                       
      #     integer i,j,k,step 
      #     REAL(kind=rc_kind) :: fac,utemp 

      for k in 1:NK
            for i in 1:NI
                  #     increase by a factor of 4 to speed spin up                        
                  #              vsouth(i,k)= 1.d-3                                       
                  vsouth[i, k] = 0e0
                  vfbcs[i, k] = vsouth[i, k] * vy[i, 1] * Jjfc[i, 0, k]
            end
      end

      #c     assign values to entrance values of s and T :arrays sent & Tent  
      for k in 1:NK
            for i in 1:NI
                  ssouth[i, k] = 33e0 - S0
                  #             ssouth[i,k]= 5.7e0 - S0                                   
                  #     ssouth is used to specify bc in sTbc_periodicew and advecn        
            end
      end
end
