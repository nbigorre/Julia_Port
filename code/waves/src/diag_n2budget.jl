function diag_n2budget(step)
   # !     =======================                                           
   # USE header
   # !     Use the Gauss' Thm form - i.e. fluxes at  top and bottom only.    
   # !     Calculate N2  (units: per s^2) and terms in EQ 6 in Thomas-Ferrari
   # !     bouyancy b in (m per s^2)                                         
   # !     vor in (per second)                                               
   # !     PV in (per second^3)                                              
   # !     zeta*b is vor*b. So zeta*b/(f H) is in (per s^2)                  
   #                                                                       
   #     implicit none 
   #     integer step,i,i2,j,k,m,n,jup,jdn,iup,idn,kup,kdn 
   #     INTEGER :: ktest(3)
   #     REAL(kind=rc_kind) :: depth,dz,dzt,dzb,dyinv,dxinv,dzinv,zbt,          &
   #    &     flocalinv,flocal,bx,by,bz,grav,adv,vol                       
   #     REAL(kind=rc_kind) :: avezf(NK),avezc(NK),avexyN2(NK),testdepth(3) 
   #     REAL(kind=rc_kind) :: b(NI,NJ,NK),Dbdt(NI,NJ,NK) 
   #     pv is in header and is common                                     
   #     REAL(kind=rc_kind) :: pvx,pvy,pvz,gravR0,dia,                                      &
   #    &     gradBxF_top,gradBxF_bot,adv_top,adv_bot,dia_top,dia_bot      
   #                                                                       
   #     WL= eps*delta*UL 
   local b = zeros(rc_kind, NI, NJ, NK)
   local Dbdt = zeros(rc_kind, NI, NJ, NK)
   local avexyN2 = zeros(rc_kind, NK)
   local avezc = zeros(rc_kind, NK)
   local avezf = zeros(rc_kind, NK)    
   local ktest = zeros(Int, 3)                                         
   local gravR0 = gpr * 10e0 / R0
   local n = 0

   #     Horizontally average                                              
   for k in 1:NK
      avexyN2[k] = 0e0
      avezf[k] = 0e0
      avezc[k] = 0e0
      for j in 1:NJ
         for i in 1:NI
            avexyN2[k] = freqN2[i, j, k] + avexyN2[k]
            avezf[k] = -zf[i, j, k] + avezf[k]
            avezc[k] = -zc[i, j, k] + avezc[k]
         end
      end
      avexyN2[k] = avexyN2[k] / rc_kind(NI * NJ)
      avezf[k] = DL * avezf[k] / rc_kind(NI * NJ)
      avezc[k] = DL * avezc[k] / rc_kind(NI * NJ)
   end
   #     avezf and avezc are positive (in m)                               

   #     Integrate vertically, down to 0.9*MLD, 0.75*MLD , and 0.4*MLD     
   local testdepth = [
      0.9e0 * mldepth,
      0.75e0 * mldepth,
      0.5e0 * mldepth
   ]

   for m in 1:3
      for k in NK:-1:1
         #     avezf and avezc are positive (refer to depth of grid level)       
         if (avezf[k-1] >= testdepth[m])
            local dzb = avezf[k-1] - testdepth[m]
            local dzt = testdepth[m] - avezf[k]
            if (dzt <= dzb)
               ktest[m] = k
            else
               ktest[m] = k - 1
            end
            break
         end
      end
   end

   #     Calculate mldN2                                                   
   #     ---------------                                                   
   for m in 1:3
      depth = 0e0
      mldN2[m] = 0e0
      for k in NK:-1:1
         if (k >= ktest[m])
            # avezf is positive (depth)  
            dz = -(avezf[k] - avezf[k-1])
            mldN2[m] = mldN2[m] + avexyN2[k] * dz
            depth = depth + dz
         end
      end
      mldN2[m] = mldN2[m] / depth

      if (step == 0)
         mldN2init[m] = mldN2[m]
      else
         mldN2[m] = mldN2[m] - mldN2init[m]
      end
   end

   #     FRONTAL TERM                                                      
   #     =============================                                     
   #     Calculate zeta*b top and bot divided by fH in per s^2     
   evalrho(rho, 0)
   vort(0)
   #@ccall "./PSOM_LIB.so".vort_(Ref(0)::Ref{Int})::Cvoid

   #     Calculate b  in m/s^2                                             
   #     -----------------------                                           
   for k in 1:NK
      for j in 1:NJ
         for i in 1:NI
            b[i, j, k] = -gravR0 * (rho[i, j, k] - R0)
         end
      end
   end

   for m in 1:3
      zbbot[m] = 0e0
   end
   zbt = 0e0
   for j in 1:NJ
      for i in 1:NI
         local flocalinv = 1e0 / (ffc[i, j] * FPAR)
         zbt = zbt - vor[i, j, NK] * b[i, j, NK] * flocalinv
         for m in 1:3
            zbbot[m] = zbbot[m] - vor[i, j, ktest[m]] * b[i, j, ktest[m]] * flocalinv
         end
      end
   end
   for m in 1:3
      zbtop[m] = zbt / rc_kind(NI * NJ) / testdepth[m]
      zbbot[m] = zbbot[m] / rc_kind(NI * NJ) / testdepth[m]
   end

   if (step == 0)
      for m in 1:3
         zbtopinit[m] = zbtop[m]
         zbbotinit[m] = zbbot[m]
      end
   else
      for m in 1:3
         zbtop[m] = zbtop[m] - zbtopinit[m]
         zbbot[m] = zbbot[m] - zbbotinit[m]
      end
   end
   #     =======================================================           
   #     ADVEC term  u \cdot  grad (q)                                     
   #     Calculate Time_integ [ volume_integ (Div ( u q )) ]   in per s^2  
   #     i.e.  Time_integ ( Volume Integ (u \cdot grad(q) ))               
   diag_pv(0)
   #     PV is in per s^3, w is in m/s, testdepth in m, f in per s         


   for k in 1:NK
      for j in 1:NJ
         for i in 0:NI+1
            if (i == 0)
               i2 = NI
            elseif (i == NI + 1)
               i2 = 1
            else
               i2 = i
            end
            pv[i, j, k] = pv1[i2, j, k] + pv2[i2, j, k] + pv3[i2, j, k]
         end
      end
   end

   for k in 1:NK
      for j in 1:NJ
         for i in 1:NI
            Dbdt[i, j, k] = -gravR0 * ((rho[i, j, k] - rhoprev[i, j, k]) / dtf + rhoadv[i, j, k]) / TL
            rhoprev[i, j, k] = rho[i, j, k]
         end
      end
   end
   #     Integrate up to ktest[m] of deepest level, i.e. m=3               

   adv = 0e0
   dia = 0e0
   for m in 1:3
      friction[m] = 0e0
      diabatic[m] = 0e0
      advecpv[m] = 0e0
      frictop[m] = 0e0
      diatop[m] = 0e0
      fricbot[m] = 0e0
      diabot[m] = 0e0
   end

   #     Integrate upto the deepest testdepth, i.e. testdepth(1)           
   for j in 1:NJ
      local jup = j + 1
      local jdn = j - 1
      local dyinv = 1e0 / (2e0 * dyM[j])
      if (j == NJ)
         jup = j
         dyinv = 1e0 / dyM[j]
      elseif (j == 1)
         jdn = j
         dyinv = 1e0 / dyM[j]
      end
      for i in 1:NI

         iup = i + 1
         idn = i - 1
         dxinv = 1e0 / (2e0 * dx)
         #     Periodic-ew boundaries                                            
         if (i == NI)
            iup = 1
         elseif (i == 1)
            idn = NI
         end

         flocalinv = 1e0 / (ffc[i, j] * FPAR)
         #     TOP                                                               
         local k = NK
         #     ADVEC TERM                                                        
         local adv_top = 0e0

         #     FRIC TERM                                                         
         #     CALCULATE bx,by,bz                                                
         bx = dxinv * (b[iup, j, k] - b[idn, j, k])
         by = dyinv * (b[i, jup, k] - b[i, jdn, k])




         gradBxF_top = bx * fricv[i, j, k] - by * fricu[i, j, k]

         dia_top = -(vor3[i, j, k] + ffc[i, j] * FPAR) * Dbdt[i, j, k]


         for m in 1:3
            k = ktest[m]

            adv_bot = -flocalinv * WL * w[i, j, k, 0] * pv[i, j, k]

            #     CALCULATE bx,by,bz                                                
            bx = dxinv * (b[iup, j, k] - b[idn, j, k])
            by = dyinv * (b[i, jup, k] - b[i, jdn, k])
            gradBxF_bot = bx * fricv[i, j, k] - by * fricu[i, j, k]

            dia_bot = -(vor3[i, j, k] + ffc[i, j] * FPAR) * Dbdt[i, j, k]
            advecpv[m] = advecpv[m] + dtf * TL * (adv_top - adv_bot) / testdepth[m]
            friction[m] = friction[m] + dtf * TL * flocalinv * (gradBxF_top - gradBxF_bot) / testdepth[m]
            frictop[m] = frictop[m] + dtf * TL * flocalinv * gradBxF_top / testdepth[m]
            fricbot[m] = fricbot[m] + dtf * TL * flocalinv * gradBxF_bot / testdepth[m]
            diabatic[m] = diabatic[m] + dtf * TL * flocalinv * (dia_top - dia_bot) / testdepth[m]
            diatop[m] = diatop[m] + dtf * TL * flocalinv * dia_top / testdepth[m]
            diabot[m] = diabot[m] + dtf * TL * flocalinv * dia_bot / testdepth[m]

         end
      end
   end

   for m in 1:3
      advecpv[m] = advecpv[m] / rc_kind(NI * NJ)
      friction[m] = friction[m] / rc_kind(NI * NJ)
      diabatic[m] = diabatic[m] / rc_kind(NI * NJ)
      frictop[m] = frictop[m] / rc_kind(NI * NJ)
      diatop[m] = diatop[m] / rc_kind(NI * NJ)
      fricbot[m] = fricbot[m] / rc_kind(NI * NJ)
      diabot[m] = diabot[m] / rc_kind(NI * NJ)
   end
end
