function write_cdf_2D_x(islice, counter_2d, n)

   #  !     ------------------------------------------------------------      
   #  USE header,ONLY : NI,NJ,NK,ntr,s,T,rho,u,v,w,consump,Tr,shear,vor,freqN2,yc,zc,fricb,stress_top_x,DL,time_seconds,dirout,rc_kind, pvt, rp, stress_top_y, pv3
   #  !     reads in a netcdf file and interpolates the data from the sigma gr
   #  !     onto a z level                                                    
   ##include "netcdf.inc"                                                                        
   #!      INCLUDE '/usr/local/include/netcdf.inc' 
   #!      INCLUDE '/sw/lib/netcdf-gfortran/include/netcdf.inc'             
   #                                                                        
   #      integer NI2,NJ2,NK2,i,j,k,n,islice,it,itr 
   #      parameter ( NI2=NI+2,NJ2=NJ+2,NK2=NK+2) 
   #      REAL(kind=rc_kind) :: trbar,fbmean 
   #      REAL(kind=rc_kind) :: svert(NK),Tvert(NK),rvert(NK),uvert(NK),         &
   #     &     vvert(NK),wvert(NK),zvert(NK),seval                          
   #                                                                        
   #      REAL(kind=4) ::  yslice(NJ),zslice(NJ,NK),conslice(ntr,NJ,NK),              &
   #           sslice(NJ,NK),Tslice(NJ,NK),rslice(NJ,NK),                   &
   #           Trslice(ntr,NJ,NK),shearslice(NJ,NK),                        &
   #           uslice (NJ,NK),vslice (NJ,NK),wslice(NJ,NK),vorslice(NJ,NK), &
   #	   ugslice(NJ,NK),vgslice(NJ,NK),                               &
   #           n2slice(NJ,NK),n2barslice(NJ,NK),zcave(NK),fricbsl(NJ,NK),   &
   #           timday,stressxsp(NJ), stressysp(NJ), pvtslice(NJ,NK), pv3slice(NJ,NK), pslice(NJ,NK)
   #!     Zonally averaged psi and by,bz                                    
   #                                                                        
   #      REAL(kind=rc_kind) ::   psiw(NJ,NK),psiv(NJ,NK),bybar(NJ,NK),          &
   #     &     bzbar(NJ,NK),vbbar(NJ,NK),wbbar(NJ,NK),psieddy(NJ,NK),       &
   #     &     rhobar(NJ,NK),n2bar(NJ,NK),drhobardz,drhodz,                 &
   #     &     wpvbar(NJ,NK),pvbar(NJ,NK),by(NJ,NK)                      
   #      REAL(kind=rc_kind) :: Feddydiv(NJ,NK),Freyndiv(NJ,NK),                 &
   #     &     cybar(ntr,NJ,NK),czbar(ntr,NJ,NK),                           &
   #     &     vcbar(ntr,NJ,NK),wcbar(ntr,NJ,NK),csfeddy(ntr,NJ,NK)         
   #                                                                        
   #      REAL(kind=4) ::  byslice(NJ,NK),bysection(NJ,NK),bzslice(NJ,NK),            &
   #     &      psivslice(NJ,NK),            &
   #     &     psiwslice(NJ,NK),vbbarsl(NJ,NK),wbbarsl(NJ,NK),              &
   #     &     psieddysl(NJ,NK),rhobarslice(NJ,NK),                         &
   #     &     wpvbarsl(NJ,NK),pvbarsl(NJ,NK),trbarsl(ntr,NJ,NK)            
   #                                                                        
   #      REAL(kind=rc_kind) ::  Feddydivsl(NJ,NK),Freyndivsl(NJ,NK),cybarsl(ntr,NJ,NK),    &
   #     &     czbarsl(ntr,NJ,NK),vcbarsl(ntr,NJ,NK),                       &
   #     &     wcbarsl(ntr,NJ,NK),csfeddysl(ntr,NJ,NK)                      
   #                                                                       
   #
   #      REAL(kind=rc_kind), dimension(0:NI+1,0:NJ+1,0:NK+1) :: ugeo, vgeo
   #      
   #      character(LEN=150) outname 
   #     
   #       INTEGER :: counter_2d
   #                                                                   
   #      integer start(3),count(3),dims(3),start2d(2),dims4(4),start4(4),  &
   #           count4(4),dimstim,counttim,count1d,dims1d,dimswind(2),       &
   #           countwind(2),startwind(2)
   #
   # integer :: iddatfile, idigit, idudx, idudy, idudz, idvbysection, idvby, idvbz, idvcon, idvcy, idvbx
   #  integer :: idvcz, idvc, idvdivfeddy, idvdivfreyn, idvdx, idvdz, idvd, idvfb, idvh, idvn2bar
   #  integer :: idvn2, idvnsq100m, idvnsq30m, idvpe,idvpsiv,idvpsiw,idvpv, idvp,idvrhbar,idvrho,idvrnk
   #  integer :: idvstrain,idvstressx,idvstressy,idvstr,idvs,idvtbar,idvtemp,idvtim,idvtr,idvt,idvu,idvvb,idvvc,idvvor,idvv, idvug, idvvg
   #  integer :: idvwb,idvwc,idvwpv,idvw,idvy,idvzsave,idvz,idwdx,idwdy,idwdz,iimday,ipos
   #  integer :: idvzave,idvshear, idvpvt, idvpv3
   #  REAL(kind=rc_kind) ::  rcode

   local yslice = zeros(Float32, NJ)
   local zslice = zeros(Float32, NJ, NK)
   local conslice = zeros(Float32, ntr, NJ, NK)
   local sslice = zeros(Float32, NJ, NK)
   local Tslice = zeros(Float32, NJ, NK)
   local rslice = zeros(Float32, NJ, NK)
   local Trslice = zeros(Float32, ntr, NJ, NK)
   local shearslice = zeros(Float32, NJ, NK)
   local uslice = zeros(Float32, NJ, NK)
   local vslice = zeros(Float32, NJ, NK)
   local wslice = zeros(Float32, NJ, NK)
   local vorslice = zeros(Float32, NJ, NK)
   local ugslice = zeros(Float32, NJ, NK)
   local vgslice = zeros(Float32, NJ, NK)
   local n2slice = zeros(Float32, NJ, NK)
   local n2barslice = zeros(Float32, NJ, NK)
   local zcave = zeros(Float32, NK)
   local fricbsl = zeros(Float32, NJ, NK)
   local stressxsp = zeros(Float32, NJ)
   local stressysp = zeros(Float32, NJ)
   local pvtslice = zeros(Float32, NJ, NK)
   local pv3slice = zeros(Float32, NJ, NK)
   local pslice = zeros(Float32, NJ, NK)
   local byslice = zeros(Float32, NJ, NK)
   local bysection = zeros(Float32, NJ, NK)
   local bzslice = zeros(Float32, NJ, NK)
   local psivslice = zeros(Float32, NJ, NK)
   local psiwslice = zeros(Float32, NJ, NK)
   local vbbarsl = zeros(Float32, NJ, NK)
   local wbbarsl = zeros(Float32, NJ, NK)
   local psieddysl = zeros(Float32, NJ, NK)
   local rhobarslice = zeros(Float32, NJ, NK)
   local wpvbarsl = zeros(Float32, NJ, NK)
   local pvbarsl = zeros(Float32, NJ, NK)
   local trbarsl = zeros(Float32, ntr, NJ, NK)



   local timday = time_seconds / 86400e0
   local rhobar = zeros(rc_kind, NJ, NK)
   local n2bar = zeros(rc_kind, NJ, NK)
   local psiv = zeros(rc_kind, NJ, NK)
   local psiw = zeros(rc_kind, NJ, NK)
   local by = zeros(rc_kind, NJ, NK)
   local bybar = zeros(rc_kind, NJ, NK)
   local bzbar = zeros(rc_kind, NJ, NK)
   local vbbar = zeros(rc_kind, NJ, NK)
   local wbbar = zeros(rc_kind, NJ, NK)
   local psieddy = zeros(rc_kind, NJ, NK)
   local wpvbar = zeros(rc_kind, NJ, NK)
   local pvbar = zeros(rc_kind, NJ, NK)
   local Feddydiv = zeros(rc_kind, NJ, NK)
   local Freyndiv = zeros(rc_kind, NJ, NK)
   local cybar = zeros(rc_kind, ntr, NJ, NK)
   local czbar = zeros(rc_kind, ntr, NJ, NK)
   local vcbar = zeros(rc_kind, ntr, NJ, NK)
   local wcbar = zeros(rc_kind, ntr, NJ, NK)



   local ugeo = OffsetArrays.zeros(rc_kind, 0:NI+1, 0:NJ+1, 0:NK+1)
   local vgeo = OffsetArrays.zeros(rc_kind, 0:NI+1, 0:NJ+1, 0:NK+1)

   diag_streamfunction(rhobar,n2bar,psiv,psiw,by,bybar,bzbar,vbbar,wbbar,psieddy,wpvbar,pvbar,Feddydiv,Freyndiv,cybar,czbar,vcbar,wcbar) 
   
   #@ccall "./PSOM_LIB.so".diag_streamfunction_(pointer(rhobar)::Ptr{rc_kind}, pointer(n2bar)::Ptr{rc_kind}, pointer(psiv)::Ptr{rc_kind}, pointer(psiw)::Ptr{rc_kind}, pointer(by)::Ptr{rc_kind}, pointer(bybar)::Ptr{rc_kind}, pointer(bzbar)::Ptr{rc_kind},
   #   pointer(vbbar)::Ptr{rc_kind}, pointer(wbbar)::Ptr{rc_kind}, pointer(psieddy)::Ptr{rc_kind}, pointer(wpvbar)::Ptr{rc_kind},
   #   pointer(pvbar)::Ptr{rc_kind}, pointer(Feddydiv)::Ptr{rc_kind}, pointer(Freyndiv)::Ptr{rc_kind}, pointer(cybar)::Ptr{rc_kind},
   #   pointer(czbar)::Ptr{rc_kind}, pointer(vcbar)::Ptr{rc_kind}, pointer(wcbar)::Ptr{rc_kind})::Cvoid


   uv_geostrophic(ugeo, vgeo)
   #@ccall "./PSOM_LIB.so".uv_geostrophic_(pointer(ugeo)::Ptr{rc_kind}, pointer(vgeo)::Ptr{rc_kind})::Cvoid

   local fricbsl = zeros(rc_kind, NJ, NK)
   for k in 1:NK
      for j in 1:NJ
         local fbmean = 0e0
         for i in 1:NI
            fbmean = fricb[i, j, k] + fbmean
         end
         fbmean = fbmean / rc_kind(NI)
         fricbsl[j, k] = fbmean
      end
   end

   for k in 1:NK
      for j in 1:NJ
         #           rhobarslice[j,k]= rhobar[j,k] 
         #             psivslice[j,k]=   psiv[j,k] 
         #             psiwslice[j,k]=   psiw[j,k] 
         #               byslice[j,k]=  bybar[j,k] 
         bysection[j, k] = by[j, k]
         #           bzslice[j,k]=  bzbar[j,k] 
         #           vbbarsl[j,k]=  vbbar[j,k] 
         #           wbbarsl[j,k]=  wbbar[j,k] 
         #         psieddysl[j,k]=psieddy[j,k] 
         #          wpvbarsl[j,k]= wpvbar[j,k] 
         #           pvbarsl[j,k]=  pvbar[j,k] 

         #            Feddydivsl[j,k]= Feddydiv[j,k] 
         #            Freyndivsl[j,k]= Freyndiv[j,k] 
         #            for it in 1:ntr 
         #               cybarsl[it,j,k]= cybar[it,j,k] 
         #               czbarsl[it,j,k]= czbar[it,j,k] 
         #               vcbarsl[it,j,k]= vcbar[it,j,k] 
         #               wcbarsl[it,j,k]= wcbar[it,j,k] 
         #               csfeddysl[it,j,k]= csfeddy[it,j,k]                      
         #            end do 

      end
   end

   local trbarsl = zeros(rc_kind, ntr, NJ, NK)
   for k in 1:NK
      for j in 1:NJ
         for itr in 1:ntr
            local trbar = 0e0
            for i in 1:NI
               trbar = trbar + Tr[itr, i, j, k, n]
            end
            trbar = trbar / rc_kind(NI)
            trbarsl[itr, j, k] = trbar
         end
      end
   end

   local i = islice
   for k in 1:NK
      zcave[k] = 0e0
      for j in 1:NJ
         zcave[k] = zcave[k] + zc[i, j, k] * DL
         zslice[j, k] = zc[i, j, k] * DL
         yslice[j] = yc[j]

         sslice[j, k] = s[i, j, k, n]
         Tslice[j, k] = T[i, j, k, n]
         rslice[j, k] = rho[i, j, k]
         uslice[j, k] = u[i, j, k, n]
         vslice[j, k] = v[i, j, k, n]
         ugslice[j, k] = ugeo[i, j, k]
         vgslice[j, k] = vgeo[i, j, k]
         #            shearslice[j,k]=shear[i,j,k]
         wslice[j, k] = w[i, j, k, n]
         vorslice[j, k] = vor[i, j, k]
         pvtslice[j, k] = pvt[i, j, k]
         pv3slice[j, k] = pv3[i, j, k]
         pslice[j, k] = rp[i, j, k]
         #            for itr in 1:ntr 
         #               conslice[itr,j,k]= consump[i,j,k,itr]
         #                Trslice[itr,j,k]=      Tr[itr,i,j,k,n]
         #            end
         #            if (k == 1)                                           
         #               dz = (zc[i,j,k+1]-zc[i,j,k])*DL                         
         #               drhodz= (rho[i,j,k+1]-rho[i,j,k])/dz                    
         #               drhobardz= (rhobar[j,k+1]-rhobar[j,k])/dz               
         #            elseif (k == NK)                   
         #               dz = (zc[i,j,k]-zc[i,j,k-1])*DL                         
         #               drhodz= (rho[i,j,k]-rho[i,j,k-1])/dz                    
         #               drhobardz= (rhobar[j,k]-rhobar[j,k-1])/dz               
         #            else                                                       
         #               dz = (zc[i,j,k+1]-zc[i,j,k-1])*DL                       
         #               drhodz= (rho[i,j,k+1]-rho[i,j,k-1])/dz                  
         #               drhobardz= (rhobar[j,k+1]-rhobar[j,k-1])/dz             
         #            end                                         
         #            n2slice[j,k] = -drhodz*9.81/R0                             
         #            n2barslice[j,k] = -drhobardz*9.81/R0                       
         n2slice[j, k] = freqN2[i, j, k]
         n2barslice[j, k] = n2bar[j, k]
      end
      zcave[k] = zcave[k] / rc_kind(NJ)
   end
   for j in 1:NJ
      stressxsp[j] = stress_top_x[div(NI, 2), j]
      stressysp[j] = stress_top_y[div(NI, 2), j]
   end


   #     write to a netcdf file  
   local outname = string("xslice_", lpad(islice, 3, "0"), ".cdf")

   if (counter_2d == 1)
      local dims = [
         NcDim("y", NJ),
         NcDim("z", NK),
         NcDim("time", 0, unlimited=true)
      ]
      local dims1d = dims[2]

      local dimswind = [dims[1], dims[3]]

      local dims4 = [NcDim("ntr", ntr), dims...]

      local dimstim = dims[3]

      local varlist = [
         NcVar("yc", dims[1], t=Float32),
         #NcVar("zcave", dims1d, t=Float32),
         NcVar("zc", dims, t=Float32),
         NcVar("day", dimstim, t=Float32),
         NcVar("stressx", dimswind, t=Float32),
         NcVar("stressy", dimswind, t=Float32),
         #NcVar("consump", dims4, t=Float32),
         #NcVar("tr", dims4, t=Float32),
         #NcVar("trbar", dims4, t=Float32),
         NcVar("s", dims, t=Float32),
         NcVar("temp", dims, t=Float32),
         NcVar("rho", dims, t=Float32),
         NcVar("u", dims, t=Float32),
         NcVar("ug", dims, t=Float32),
         NcVar("v", dims, t=Float32),
         NcVar("vg", dims, t=Float32),
         NcVar("w", dims, t=Float32),
         NcVar("vor", dims, t=Float32),
         #NcVar("shear", dims, t=Float32),
         NcVar("n2", dims, t=Float32),
         #NcVar("fricb", dims, t=Float32),
         NcVar("pvt", dims, t=Float32),
         NcVar("pv3", dims, t=Float32),
         NcVar("pre", dims, t=Float32),
         #NcVar("rhobar", dims, t=Float32),
         #NcVar("n2bar", dims, t=Float32),
         #NcVar("psivbar", dims, t=Float32),
         #NcVar("psiwbar", dims, t=Float32),
         NcVar("by", dims, t=Float32),
         #NcVar("bybar", dims, t=Float32),
         #NcVar("bzbar", dims, t=Float32),
         #NcVar("vbbar", dims, t=Float32),
         #NcVar("wbbar", dims, t=Float32),
         #NcVar("psieddy", dims, t=Float32),
         #NcVar("pvbar", dims, t=Float32),
         #NcVar("wpvbar", dims, t=Float32),
         #NcVar("divfreyn", dims, t=Float32),
         #NcVar("divfreddy", dims, t=Float32),
         #NcVar("cybar", dims4, t=Float32),
         #NcVar("czbar", dims4, t=Float32),
         #NcVar("vcbar", dims4, t=Float32),
         #NcVar("wcbar", dims4, t=Float32),
         #NcVar("psieddytracer", dims4, t=Float32)

      ]

      NetCDF.create(string(dirout, "/", outname), varlist, mode=NC_CLASSIC_MODEL)
   end

   local counttim = [1]
   local count1d = [NK]

   local count = [NJ, NK, 1]
   local countwind = [NJ, 1]
   local count4 = [ntr, NJ, NK, 1]
   local start2d = [1, 1]
   local start = [1, 1, counter_2d]
   local startwind = [1, counter_2d]
   local start4 = [1, 1, 1, counter_2d]


   local ncfile = NetCDF.open(string(dirout, "/", outname), mode=NC_WRITE)
   if (counter_2d == 1)
      NetCDF.putvar(ncfile["yc"], yslice, start=[start[1]], count=[count[1]])
      #NetCDF.putvar(ncfile["zcave"], zcave, start= start(1), count=count1d) 
      NetCDF.putvar(ncfile["zc"], zslice, start=start, count=count)
   end
   NetCDF.putvar(ncfile["day"], [timday], start=[start[3]], count=counttim)
   NetCDF.putvar(ncfile["stressx"], stressxsp, start=startwind, count=countwind)
   NetCDF.putvar(ncfile["stressy"], stressysp, start=startwind, count=countwind)
   #NetCDF.putvar(ncfile["consump"],  conslice, start= start4, count= count4) 
   #NetCDF.putvar(ncfile["tr"],  Trslice, start= start4, count= count4) 
   #NetCDF.putvar(ncfile["trbar"],  Trbarsl, start= start4, count= count4) 
   NetCDF.putvar(ncfile["s"], sslice, start=start, count=count)
   NetCDF.putvar(ncfile["temp"], Tslice, start=start, count=count)
   NetCDF.putvar(ncfile["rho"], rslice, start=start, count=count)
   NetCDF.putvar(ncfile["u"], uslice, start=start, count=count)
   NetCDF.putvar(ncfile["ug"], ugslice, start=start, count=count)
   NetCDF.putvar(ncfile["v"], vslice, start=start, count=count)
   NetCDF.putvar(ncfile["vg"], vgslice, start=start, count=count)
   NetCDF.putvar(ncfile["w"], wslice, start=start, count=count)
   NetCDF.putvar(ncfile["vor"], vorslice, start=start, count=count)
   #NetCDF.putvar(ncfile["shear"],  shearslice, start= start, count= count) 
   NetCDF.putvar(ncfile["n2"], n2slice, start=start, count=count)
   #NetCDF.putvar(ncfile["fricb"],  fricbsl, start= start, count= count) 
   NetCDF.putvar(ncfile["pvt"], pvtslice, start=start, count=count)
   NetCDF.putvar(ncfile["pv3"], pv3slice, start=start, count=count)
   NetCDF.putvar(ncfile["pre"], pslice, start=start, count=count)

   #NetCDF.putvar(ncfile["rhobar"],  rhobarslice, start= start, count= count) 
   #NetCDF.putvar(ncfile["n2bar"],  n2barslice, start= start, count= count) 
   #NetCDF.putvar(ncfile["psivbar"],  psivslice, start= start, count= count) 
   #NetCDF.putvar(ncfile["psiwbar"],  psiwslice, start= start, count= count) 
   NetCDF.putvar(ncfile["by"], bysection, start=start, count=count)
   #NetCDF.putvar(ncfile["bybar"],  byslice, start= start, count= count) 
   #NetCDF.putvar(ncfile["bzbar"],  bzslice, start= start, count= count) 
   #NetCDF.putvar(ncfile["vbbar"],  vbbarsl, start= start, count= count) 
   #NetCDF.putvar(ncfile["wbbar"],  wbbarsl, start= start, count= count) 
   #NetCDF.putvar(ncfile["psieddy"],  psieddysl, start= start, count= count) 
   #NetCDF.putvar(ncfile["pvbar"],  pvbarsl, start= start, count= count) 
   #NetCDF.putvar(ncfile["wpvbar"],  wpvbarsl, start= start, count= count) 

   #NetCDF.putvar(ncfile["divfreyn"],  Freyndivsl, start= start, count= count) 
   #NetCDF.putvar(ncfile["divfeddy"],  Feddydivsl, start= start, count= count) 

   #NetCDF.putvar(ncfile["cybar"],  cybarsl, start=start4, count=count4) 
   #NetCDF.putvar(ncfile["czbar"],  czbarsl, start=start4, count=count4) 
   #NetCDF.putvar(ncfile["vcbar"],  vcbarsl, start=start4, count=count4) 
   #NetCDF.putvar(ncfile["wcbar"],  wcbarsl, start=start4, count=count4) 
   #NetCDF.putvar(ncfile["psieddytracer"],  csfeddysl, start=start4, count=count4)   

   NetCDF.sync(ncfile)
end
