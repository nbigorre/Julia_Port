function diag_streamfunction(rhobar, n2bar, psiv, psiw, by, bybar, bzbar, vbbar, wbbar, psieddy, wpvbar, pvbar, Feddydiv, Freyndiv, cybar, czbar, vcbar, wcbar)
   #     --------------                                                    
   #  USE header
   #!     Solve psi_yy + psi_zz = vor in 2-D                                
   #!     BC psi =0 on all boundaries                                       
   #                                                                        
   #      implicit NONE 
   #      integer i,j,k,n,jup,jdn,kup,kdn
   #      REAL(kind=rc_kind) :: fcor(NJ),wzbar(NJ,NK),vybar(NJ),vbar(NJ,NK),     &
   #     &     wbar(NJ,NK),psiw(NJ,NK),psiv(NJ,NK),bybar(NJ,NK), &
   #     &     bzbar(NJ,NK),bbar(NJ,NK),zcbar(NJ,NK),ubar(NJ,NK),           &
   #     &     b(NI,NJ,NK),by(NJ,NK),bz(NI,NJ,NK),                       &
   #     &     vbbar(NJ,NK),wbbar(NJ,NK),psieddy(NJ,NK),rhobar(NJ,NK),      &
   #     &     wpvbar(NJ,NK),pvbar(NJ,NK),                    &
   #     &     pv1bar(NJ,NK),pv2bar(NJ,NK),     &
   #     &     pv3bar(NJ,NK),n2bar(NJ,NK)                                   
   #      REAL(kind=rc_kind) :: Feddydiv(NJ,NK),Freyndiv(NJ,NK),                 &
   #     &     cybar(ntr,NJ,NK),czbar(ntr,NJ,NK),                           &
   #     &     cbar(ntr,NJ,NK),                                             &
   #     &     vcbar(ntr,NJ,NK),wcbar(ntr,NJ,NK),csfeddy(ntr,NJ,NK)         
   #      REAL(kind=rc_kind) :: psivprev,psiwprev,dyc,dz,vp,wp,bp,cp,par,        &
   #     &     rhomean,dyinv                                                
   #           
   local fcor = zeros(rc_kind, NJ)
   local wzbar = zeros(rc_kind, NJ, NK)
   local vybar = zeros(rc_kind, NJ)
   local vbar = zeros(rc_kind, NJ, NK)
   local wbar = zeros(rc_kind, NJ, NK)
   local bbar = zeros(rc_kind, NJ, NK)
   local zcbar = zeros(rc_kind, NJ, NK)
   local ubar = zeros(rc_kind, NJ, NK)
   local b = zeros(rc_kind, NI, NJ, NK)
   local pv1bar = zeros(rc_kind, NJ, NK)
   local pv2bar = zeros(rc_kind, NJ, NK)
   local pv3bar = zeros(rc_kind, NJ, NK)
   local cbar = zeros(rc_kind, ntr, NJ, NK)
   local n = 0

   evalrho(rho, 0)
   #@ccall "./PSOM_LIB.so".evalrho_(pointer(rho)::Ptr{rc_kind}, Ref(0)::Ref{Int})::Cvoid
   diag_pv(n)
   #@ccall "./PSOM_LIB.so".diag_pv_(Ref(n)::Ref{Int})::Cvoid

   local rhomean = 0e0
   for k in 1:NK
      for j in 1:NJ
         wzbar[j, k] = 0e0
         for i in 1:NI
            wzbar[j, k] = wzbar[j, k] + wz[i, j, k]
            rhomean = rho[i, j, k] + rhomean
         end
         wzbar[j, k] = wzbar[j, k] / rc_kind(NI)
      end
   end
   rhomean = rhomean / rc_kind(NI * NJ * NK)

   #     Calculate b  in m/s^2                                             
   for k in 1:NK
      for j in 1:NJ
         for i in 1:NI
            #               b[i,j,k]= -gpr*10e0*(rho[i,j,k]-rhomean)/rhomean      
            b[i, j, k] = -gpr * 10e0 * (rho[i, j, k] - R0) / R0
         end
      end
   end
   #      calculate by on section i=NI/2
   local i = div(NI, 2)
   for k in 1:NK
      for j in 1:NJ
         local jup = j + 1
         local jdn = j - 1
         if (j == 1)
            jdn = j
         end
         if (j == NJ)
            jup = j
         end
         local dyc = (yc[jup] - yc[jdn]) * 1e3
         by[j, k] = (b[i, jup, k] - b[i, jdn, k]) / dyc
      end
   end

   #     Average zonally                                                   
   #     =================                                                 
   for k in 1:NK
      for j in 1:NJ

         ubar[j, k] = 0e0
         vbar[j, k] = 0e0
         wbar[j, k] = 0e0
         bbar[j, k] = 0e0
         rhobar[j, k] = 0e0
         n2bar[j, k] = 0e0
         zcbar[j, k] = 0e0
         for it in 1:ntr
            cbar[it, j, k] = 0e0
         end

         pv1bar[j, k] = 0e0
         pv2bar[j, k] = 0e0
         pv3bar[j, k] = 0e0
         pvbar[j, k] = 0e0

         for i in 1:NI
            ubar[j, k] = ubar[j, k] + u[i, j, k, n] * UL
            vbar[j, k] = vbar[j, k] + v[i, j, k, n] * UL
            wbar[j, k] = wbar[j, k] + w[i, j, k, n] * WL
            bbar[j, k] = bbar[j, k] + b[i, j, k]
            rhobar[j, k] = rhobar[j, k] + rho[i, j, k]
            n2bar[j, k] = n2bar[j, k] + freqN2[i, j, k]
            zcbar[j, k] = zcbar[j, k] + zc[i, j, k]
            for it in 1:ntr
               cbar[it, j, k] = cbar[it, j, k] + Tr[it, i, j, k, n]
            end

            pv1bar[j, k] = pv1bar[j, k] + pv1[i, j, k]
            pv2bar[j, k] = pv2bar[j, k] + pv2[i, j, k]
            pv3bar[j, k] = pv3bar[j, k] + pv3[i, j, k]
            pvbar[j, k] = pv1bar[j, k] + pv2bar[j, k] + pv3bar[j, k]
         end
         ubar[j, k] = ubar[j, k] / rc_kind(NI)
         vbar[j, k] = vbar[j, k] / rc_kind(NI)
         wbar[j, k] = wbar[j, k] / rc_kind(NI)
         bbar[j, k] = bbar[j, k] / rc_kind(NI)
         rhobar[j, k] = rhobar[j, k] / rc_kind(NI)
         n2bar[j, k] = n2bar[j, k] / rc_kind(NI)
         zcbar[j, k] = DL * zcbar[j, k] / rc_kind(NI)
         for it in 1:ntr
            cbar[it, j, k] = cbar[it, j, k] / rc_kind(NI)
         end

         pv1bar[j, k] = pv1bar[j, k] / rc_kind(NI)
         pv2bar[j, k] = pv2bar[j, k] / rc_kind(NI)
         pv3bar[j, k] = pv3bar[j, k] / rc_kind(NI)
         pvbar[j, k] = pvbar[j, k] / rc_kind(NI)
      end
   end

   for j in 1:NJ
      vybar[j] = 0e0
      for i in 1:NI
         vybar[j] = vybar[j] + vy[i, j]
      end
      vybar[j] = vybar[j] / rc_kind(NI)
      fcor[j] = FPAR * ffc[div(NI, 2), j]
   end

   for k in 1:NK
      for j in 1:NJ
         vbbar[j, k] = 0e0
         wbbar[j, k] = 0e0
         wpvbar[j, k] = 0e0
         for it in 1:ntr
            vcbar[it, j, k] = 0e0
            wcbar[it, j, k] = 0e0
         end
         for i in 1:NI
            vp = UL * v[i, j, k, n] - vbar[j, k]
            wp = WL * w[i, j, k, n] - wbar[j, k]
            bp = b[i, j, k] - bbar[j, k]
            for it in 1:ntr
               cp = Tr[it, i, j, k, n] - cbar[it, j, k]
               vcbar[it, j, k] = vcbar[it, j, k] + vp * cp
               wcbar[it, j, k] = wcbar[it, j, k] + wp * cp
            end

            vbbar[j, k] = vbbar[j, k] + vp * bp
            wbbar[j, k] = wbbar[j, k] + wp * bp

            wpvbar[j, k] = wpvbar[j, k] + w[i, j, k, n] * WL * (pv1[i, j, k] + pv2[i, j, k] + pv3[i, j, k])
         end
         vbbar[j, k] = vbbar[j, k] / rc_kind(NI)
         wbbar[j, k] = wbbar[j, k] / rc_kind(NI)
         wpvbar[j, k] = wpvbar[j, k] / rc_kind(NI)
         for it in 1:ntr
            vcbar[it, j, k] = vcbar[it, j, k] / rc_kind(NI)
            wcbar[it, j, k] = wcbar[it, j, k] / rc_kind(NI)
         end
      end
   end
   #     Y-derivative                                                      
   #     ---------------                                                   
   for k in 1:NK
      for j in 1:NJ
         local jup = j + 1
         local jdn = j - 1
         if (j == 1)
            jdn = j
         end
         if (j == NJ)
            jup = j
         end
         local dyc = (yc[jup] - yc[jdn]) * 1e3
         bybar[j, k] = (bbar[jup, k] - bbar[jdn, k]) / dyc
         #     b_y is in per s^2                                                 
         for it in 1:ntr
            cybar[it, j, k] = (cbar[it, jup, k] - cbar[it, jdn, k]) / dyc
         end
         #     c_y is in per m, if c is non-dim                                  
      end
   end

   #     Z-derivative                                                      
   #     -----------------                                                 
   for k in 1:NK
      local kup = k + 1
      local kdn = k - 1
      if (k == 1)
         kdn = k
      end
      if (k == NK)
         kup = k
      end
      for j in 1:NJ
         # already multipl by DL      
         local dz = zcbar[j, kup] - zcbar[j, kdn]
         bzbar[j, k] = (bbar[j, kup] - bbar[j, kdn]) / dz
         for it in 1:ntr
            czbar[it, j, k] = (cbar[it, j, kup] - cbar[it, j, kdn]) / dz
         end
      end
   end

   local par = 1e-3
   for k in 1:NK
      for j in 1:NJ
         psieddy[j, k] = par * (par * vbbar[j, k] * bzbar[j, k] - wbbar[j, k] * bybar[j, k] / par) / (bybar[j, k] * bybar[j, k] + par * par * bzbar[j, k] * bzbar[j, k])
         #     Eddy stream function based on tracer blows up where cy,cz are zero
         #            for it in 1:ntr                                                
         #               csfeddy[it,j,k]= par*(par*vcbar[it,j,k]*czbar[it,j,k] - wcbar[it,j,k] *cybar[it,j,k]/par) / (cybar[it,j,k]*cybar[it,j,k]+par*par*czbar[it,j,k]*czbar[it,j,k])              
         #            end
      end
   end

   local n = 0
   #     WL = EPS*delta*UL                                                 
   for j in 1:NJ
      local psivprev = 0e0
      for k in 1:NK
         psiv[j, k] = psivprev - DL * vbar[j, k] / wzbar[j, k]
         psivprev = psiv[j, k]
      end
   end

   for k in 1:NK
      local psiwprev = 00
      for j in 1:NJ
         psiw[j, k] = psiwprev + LEN * wbar[j, k] / vybar[j]
         psiwprev = psiw[j, k]
      end
   end

   #     Calculate divergence of F_reynolds and F_eddy (Residual is the sum
   #     ----------------------------------------------                    
   for k in 1:NK
      local kup = k + 1
      local kdn = k - 1
      if (k == 1)
         kdn = k
      end
      if (k == NK)
         kup = k
      end
      for j in 1:NJ
         local jup = j + 1
         local jdn = j - 1
         if (j == 1)
            jdn = j
         end
         if (j == NJ)
            jup = j
         end
         local dyc = (yc[jup] - yc[jdn]) * 1e3
         # already multipl by DL      
         dz = zcbar[j, kup] - zcbar[j, kdn]
         Freyndiv[j, k] = (vbbar[jup, k] - vbbar[jdn, k]) / dyc + (wbbar[j, kup] - wbbar[j, kdn]) / dz
         Feddydiv[j, k] = (psieddy[jup, k] * bzbar[jup, k] - psieddy[jdn, k] * bzbar[jdn, k]) / dyc - (psieddy[j, kup] * bybar[j, kup] - psieddy[j, kdn] * bybar[j, kdn]) / dz
      end
   end
end

