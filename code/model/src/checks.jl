function checks()

   #  !     ------------------                                                
   #  USE header
   #  integer i,j,k 
   #  integer n1, n2, n3
   #  real r1
   #
   # !     Checks and writes a few things                                    


   for j in 1:NJ
      for i in 0:NI
         if (abs(D[i+1, j] - D[i, j]) > 1e-4)
            println(i, " ", j, " ", D(i + 1, j), " ", D(i, j))
            println("need to modify rpevalgrad (rgradients) for slope in x direcn: grpifc needs to be modified")
            println("need to modify velbc_periodicew")
            exit(1)
         end
      end
   end

   for j in 0:NJ+1
      for i in 0:NI+1
         if ((uy[i, j] != 0e0) || (vx[i, j] != 0e0))
            println("Need to modify biharmonic for curvi grid")
            exit(1)
         end
      end
   end

   if (rect == (false)) #  need to modify grpifc,grpjfc to contain cross diff terms       
      println("modify grpifc,grpjfc, stop in checks")
      exit(1)
   end


   if (verbose)
      println("                                            ")
      println("############################################")
      println("#                                           ")
      println("#----------------  MODEL SETTING     ")
      println("#                                           ")


      println("#   time step   x number of time steps  = length of simulation ")
      #println(" # ",dtime_dim," sec  x      ",nsteps,"           = ")

      local r1 = nsteps * dtf * TL / 3600e0
      if (r1 > 24)
         local n1 = div(r1,24e0)
         r1 = r1 - n1 * 24e0
         println("# ", dtime_dim, " sec  x      ", nsteps, "           = ", n1, " days and ", r1, " h")
      else
         println("# ", dtime_dim, " sec  x      ", nsteps, "           = ", nsteps * dtf * TL / 3600e0, "h")
      end

      println("#                                           ")

      @static if (cppdefs.rhoonly)

         println("# only rho is used.")
      else
         println("# s and T are used.")
      end

      println("#                                           ")

      if (fplane == 1)
         println("# f-plane at ", phi0deg, " deg ")
      else
         println("# f varies with latitude, centered on ", phi0deg, " deg")
      end

      if (fnhhy >= 0.9)
         println("# non-hydrostatic simulation ")
      else
         println("# hydrostatic approximation is used ")

      end

      if (use_shchepetkin != 0)
         println("# Shchepetkin scheme ENABLED")
      else
         if (lv_flat_bottom != 0)
            println("# Shchepetkin scheme DISABLED")
            println("# Song scheme ENABLED")
         else
            println("# Shchepetkin scheme DISABLED")
            println("------")
            println("Error: You attempt to run a simulation with a non-flat bottom with the Song scheme.              ")
            println("         To use a baroclinic pressure term computation scheme that is accurate on sloping bottom,")
            println("          use Shchepetkin scheme: in namelist, use_Shchepetkin=.TRUE.                            ")
            exit(1)
         end
      end

      println("# horizontal diffusion coefficients (m2.s-1) :", Kx, " ", Ky)

      @static if (cppdefs.fixed_bottom_thickness)
         if (bottom_linear_drag)
            println("# bottom drag is linear, with coefficient RR=", RR)
         else
            println("# bottom drag not linear: not supported.")
            exit(1)
            #println("# bottom drag is quadratic, with coefficient RR=",RR)
         end
      else
         if (abs(RR) > 1e-13 && (lv_flat_bottom == 0))
            println("------")
            println("Error: Attempt to use bottom friction on a sloping topography without keeping the height of the lowermost cell.")
            println("       You have to modify at least one of these parameters:")
            println("        - fixed_bottom_thickness in inc/cppdefs.h (currently undef)")
            println("        - RR in namelist (currently non zero)")
            println("        - lv_flat_bottom in namelist (currently .FALSE.)")
            exit(1)
         end
      end

      println("#                                            ")
      println("#----------------  NONDIMENSIONAL NUMBERS    ")
      println("#                                            ")

      println("# EPS     = ", EPS)
      println("# delta   = ", delta)
      println("# qpr     = ", qpr)
      println("# lambda  = ", lambda)
      println("# beta    = ", beta)
      println("# dtf     = ", dtf)
      println("# WL      = ", WL)
      println("# R0      = ", R0)
      println("# UL      = ", UL)
      println("# P1      = ", P1)
      println("# HL      = ", HL)
      println("# HDL     = ", HDL)


      println("#                                           ")
      println("#----------------  GRID CHARACTERISTICS     ")
      println("#                                           ")
      if (rect)
         println("# rectangular grid")
      else
         println("# non-rectangular grid !")
      end
      if (periodicew)
         println("# periodic in the EW direction")
      else
         println("# non periodic in the EW direction (not supported!)")
      end


      println("# NI, NJ, NK: ", NI, " ", NJ, " ", NK)
      println("# delta_x: ", dx, " m")
      println("# delta_y: ", dy, " m")


      println("#                                           ")

      println("# top cell thickness (m)      : ", dztop)

      @static if (cppdefs.fixed_bottom_thickness)
         println("# bottom cell thickness (m)   : ", dzbot)
      else
         println("# bottom cell height       : variable")
      end

      println("# vertical grid stretching factor   : ", pfac)

      if (lv_flat_bottom != 0)
         println("# flat topography, depth: ", total_depth, " m")
      else
         println("# sloping topography, varying between: ", minimum(zf(:, :, 0)) * DL, " m and ", maximum(zf(:, :, 0)) * DL, " m")
      end

      println("#                                           ")
      println("#-----------  PARTICLES CHARACTERISTICS     ")
      println("#                                           ")


      @static if (cppdefs.allow_particle)
         println("# particles ENABLED.")
         println("# number of particles   : ", NPR)
         println("# initial time step     : ", ini_particle_time)
         println("# frequency of output   : ", parti_outfreq)
         println("# number of output files: ", parti_file_num)
         println("# pcx, pcy, pcz, pcr: ", pcx, " ", pcy, " ", pcz, " ", pcr)

      else
         println("# particles DISABLED. ")
      end


      println("#                                           ")
      println("#----------------  OUTPUT CHARACTERISTICS     ")
      println("#                                           ")


      @static if (cppdefs.file_output)

         println("# output ENABLED, in  ", dirout)


         @static if (cppdefs.file_output_cdf)
            println("# output in cdf format ENABLED.")
         else
            println("# output in cdf format DISABLED.")
         end

         @static if (cppdefs.file_output_bin)
            println("# output in bin format ENABLED.")
         else
            println("# output in bin format DISABLED.")
         end


         println("# output 1D frequency: ", out1d_int, "steps")
         println("# output 2D frequency: ", out2d_int, "steps")
         println("# output 3D frequency: ", out3d_int, "steps")


      else

         println("# output DISABLED.")

      end

      println("#                                           ")
      println("#----------------  USER VARIABLES           ")
      println("#                                           ")
      println("# user1 : ", user1)
      println("# user2 : ", user2)
      println("# user3 : ", user3)
      println("# user4 : ", user4)
      println("# user5 : ", user5)
      println("# user6 : ", user6)
      println("# user7 : ", user7)
      println("# user8 : ", user8)


      println("#                                           ")
      println("############################################")
      println("                                            ")
      println("                                            ")



      println("                                            ")
      println("                                            ")
      println("############################################")
      println("#                                           ")
      println("#----------------  INITIAL SETTING     ")
      println("#                                           ")
      println("#   u lies between: ", minimum(u[:, :, :, 0]), " and ", maximum(u[:, :, :, 0]))
      println("#   v lies between: ", minimum(v[:, :, :, 0]), " and ", maximum(v[:, :, :, 0]))
      println("#   w lies between: ", minimum(w[:, :, :, 0]), " and ", maximum(w[:, :, :, 0]))

      @static if (cppdefs.rhoonly)
         println("# rho lies between: ", minimum(rho[:, :, :]), " and ", maximum(rho[:, :, :]), " (only rho)")

      else
         println("#   s lies between: ", minimum(s[:, :, :, 0]), " and ", maximum(s[:, :, :, 0]))
         println("#   T lies between: ", minimum(T[:, :, :, 0]), " and ", maximum(T[:, :, :, 0]))
         println("# rho lies between: ", minimum(rho[:, :, :]), " and ", maximum(rho[:, :, :]), " (rho from s and T)")
      end
      println("#   h lies between: ", minimum(h[:, :]), " and", maximum(h[:, :]))
      println("#                                           ")
      println("############################################")











   end # verobse




end
