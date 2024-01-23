using NetCDF

function write_cdf_1D_mooring(imooring, jmooring, counter_1d)

   local ugeo = OffsetArrays.zeros(rc_kind, 0:NI+1, 0:NJ+1, 0:NK+1)
   local vgeo = OffsetArrays.zeros(rc_kind, 0:NI+1, 0:NJ+1, 0:NK+1)

   local zprof = zeros(Float32, NK)
   local n2prof = zeros(Float32, NK)
   local uprof = zeros(Float32, NK)
   local vprof = zeros(Float32, NK)
   local wprof = zeros(Float32, NK)
   local ugeoprof = zeros(Float32, NK)
   local vgeoprof = zeros(Float32, NK)
   local vorprof = zeros(Float32, NK)
   local rhoprof = zeros(Float32, NK)


   local i = imooring
   local j = jmooring
   for k in 1:NK
      rhoprof[k] = rho[i, j, k]
      zprof[k] = zc[i, j, k] * DL
      uprof[k] = u[i, j, k, 0]
      vprof[k] = v[i, j, k, 0]
      wprof[k] = w[i, j, k, 0]
      ugeoprof[k] = ugeo[i, j, k]
      vgeoprof[k] = vgeo[i, j, k]
      vorprof[k] = vor[i, j, k]
   end

   for k in 1:NK
      if (k == 1)
         dz = (zprof[k+1] - zprof[k])
         drhodz = (rhoprof[k+1] - rhoprof[k]) / dz
      elseif (k == NK)
         dz = (zprof[k] - zprof[k-1])
         drhodz = (rhoprof[k] - rhoprof[k-1]) / dz
      else
         dz = (zprof[k+1] - zprof[k-1])
         drhodz = (rhoprof[k+1] - rhoprof[k-1]) / dz
      end
      n2prof[k] = -drhodz * 9.81 / (rhoprof[k] + 1000.0)
   end

   #     write to a netcdf file                                            
   local outname = string("mooring_", lpad(imooring, 3, "0"), "_", lpad(jmooring, 3, "0"), ".cdf")

   if (counter_1d == 1)
      #idDatFile =  nccre(TRIM(dirout)//outname,NCCLOB,rcode) 
      local dims = [
         NcDim("z", NK),
         NcDim("time", 0, unlimited=true)
      ]
      local varlist = [
         NcVar("zc", [dims[1]], t=Float32),
         NcVar("rho", dims, t=Float32),
         NcVar("n2", dims, t=Float32),
         NcVar("u", dims, t=Float32),
         NcVar("v", dims, t=Float32),
         NcVar("ugeo", dims, t=Float32),
         NcVar("vgeo", dims, t=Float32),
         NcVar("w", dims, t=Float32),
         NcVar("vor", dims, t=Float32),
      ]
      NetCDF.create(string(dirout, "/", outname), varlist, mode=NC_CLASSIC_MODEL)


   end
   local ncfile = NetCDF.open(string(dirout, "/", outname), mode=NC_WRITE)

   local count = [NK, 1]
   local start = [1, counter_1d]
   if (counter_1d == 1)
      NetCDF.putvar(ncfile["zc"], zprof, start=[start[1]], count=[count[1]])
   end

   NetCDF.putvar(ncfile["rho"], rhoprof, start=start, count=count)
   NetCDF.putvar(ncfile["n2"], n2prof, start=start, count=count)
   NetCDF.putvar(ncfile["u"], uprof, start=start, count=count)
   NetCDF.putvar(ncfile["v"], vprof, start=start, count=count)
   NetCDF.putvar(ncfile["ugeo"], ugeoprof, start=start, count=count)
   NetCDF.putvar(ncfile["vgeo"], vgeoprof, start=start, count=count)
   NetCDF.putvar(ncfile["w"], wprof, start=start, count=count)
   NetCDF.putvar(ncfile["vor"], vorprof, start=start, count=count)

   NetCDF.sync(ncfile)

   return
end

#=
! subroutine uv_geostrophic(ugeo,vgeo)
! !     ---------------------------------------------                     
!   !   Finds the geostrophic velocities ugeo,vgeo
!   USE header
!   USE rpgrads
! !     modified for periodicew bc                                        
! !     ---------------------------                                       
! !     Sets up the initial velocities so that they are in geostrophic bal
! !     with the initial density field.                                   
! !      implicit logical (a-z)                                           
!       integer i,j,k,n,imax,jmax,kmax,m 
!       double precision uxi,vyj,hxi,heta,hx,hy,px,py,ujfc 
!       double precision res,resmax
!       double precision ainv,be2,fac2,wzsk,wxsk,wysk,Jack,pxi,peta,      &
!      &     psig,pgrad,con,pz                                            
!       real*8 ugeo(0:NI+1,0:NJ+1,0:NK+1),vgeo(0:NI+1,0:NJ+1,0:NK+1)
!                                                                         
!       call rpevalgrad(n) 
!                                                                         
!       kaph1= 1.d0 -kappah 
!       con = 1.0 -qpr 
!                                                                         
!       do j=1,NJ 
!          do i=1,NI 
!             hxi= 0.5d0*( h(i+1,j)-h(i-1,j) ) 
!             heta= 0.5d0*( h(i,j+1)-h(i,j-1) ) 
!             hx= ux(i,j)*hxi +vx(i,j)*heta 
!             hy= uy(i,j)*hxi +vy(i,j)*heta 
!             do k=1,NK 
!                pxi= 0.5d0*(p(i+1,j,k)-p(i-1,j,k)) 
!                peta= 0.5d0*(p(i,j+1,k)-p(i,j-1,k)) 
!                psig= 0.5d0*(p(i,j,k+1)-p(i,j,k-1)) 
!                px= (ux(i,j)*pxi +vx(i,j)*peta +wx(i,j,k)*psig) 
!                py= (uy(i,j)*pxi +vy(i,j)*peta +wy(i,j,k)*psig) 
!                                                                         
!                ugeo(i,j,k) = - (qpr*py +drpy(i,j,k)+gpr*hy)/             &
!                     (ffc(i,j))                                          
!                vgeo(i,j,k) = (qpr*px +drpx(i,j,k) +gpr*hx)/              &
!                     (ffc(i,j))                                          
!             end do
!          end do 
!          do k=1,NK 
!             ugeo(0,j,k)= ugeo(NI,j,k) 
!             vgeo(0,j,k)= vgeo(NI,j,k) 
!             ugeo(NI+1,j,k) = ugeo(1,j,k) 
!             vgeo(NI+1,j,k) = vgeo(1,j,k) 
!          end do 
!       end do 

!       return
!       end
=#


