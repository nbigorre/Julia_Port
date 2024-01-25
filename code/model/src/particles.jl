#include "cppdefs.h"
@static if (cppdefs.allow_particle)

module particles
#  USE header, ONLY : NI,NJ,NK,uf,vf,wf,Jifc,Jjfc,J2d,ux,vy,NPR,wz,PI,dtf,vor,shear,rho,strain,zf,&
#                    &s,parti_file_num,DL,rc_kind, pcx, pcy, pcz, pcr, dirout
#
#  ! define the class for particles
#  TYPE particle
#     REAL(kind=rc_kind) ::  i,j,k,x,y,z,u,v,w,s,t,u0,v0,w0,id,vor,strain,shear,rho,time
#  END TYPE particle
#
#  TYPE (particle), DIMENSION(:), ALLOCATABLE :: parti
#  REAL(kind=rc_kind) ::  dz,swap1,swap2,swap3
#  INTEGER,ALLOCATABLE :: file_id(:)
#  INTEGER :: NPR_eachfile
#  CHARACTER(len=3) :: file_id_char 
#
#  PRIVATE   :: NPR_eachfile, file_id_char, dz, swap1, swap2, swap3, file_id
#  PUBLIC    :: parti
#
#CONTAINS
mutable struct Particle
   i::rc_kind
   j::rc_kind
   k::rc_kind
   x::rc_kind
   y::rc_kind
   z::rc_kind
   u::rc_kind
   v::rc_kind
   w::rc_kind
   s::rc_kind
   t::rc_kind
   u0::rc_kind
   v0::rc_kind
   w0::rc_kind
   id::rc_kind
   vor::rc_kind
   strain::rc_kind
   shear::rc_kind
   rho::rc_kind
   time::rc_kind
 end

  NPR_eachfile::Int = 0
  fileList = zeros(String, parti_file_num)
  parti = zeros(Particle, NPR)
  function open_parti_files()

    if (mod(npr,parti_file_num) != 0) THEN
      println("Error: Please make sure NPR/file_num = integer in mod_particles.f90")
      println("Stop model")
      exit(1)
    end

    global NPR_eachfile = div(NPR,parti_file_num) #the particle number in each file

    #open files
    local dims = [
         NcDim("p_num", NPR_eachfile),
         NcDim("step", 0, unlimited=true)
      ]
      local varlist = [
         NcVar("i", dims, t=Float64),
         NcVar("j", dims, t=Float64),
         NcVar("k", dims, t=Float64),
         NcVar("z", dims, t=Float64),
         NcVar("u", dims, t=Float64),
         NcVar("v", dims, t=Float64),
         NcVar("w", dims, t=Float64),
         NcVar("rho", dims, t=Float64),
         NcVar("s", dims, t=Float64),
         NcVar("t", dims, t=Float64),
         NcVar("vor", dims, t=Float64),
         NcVar("shear", dims, t=Float64),
         NcVar("strain", dims, t=Float64),
      ]
    for fi in 1:parti_file_num
      local filename = string("op.parti-",lpad(fi,3,'0'),".cdf")
      fileList[fi] = string(dirout, "/", filename)
      NetCDF.create(fileList[fi], varlist, mode=NC_CLASSIC_MODEL)
    end
   end

  function save_parti(step)

    println("SAVE PARTICLES")

       
       #save limited variables
       local start = [1, step]
       local count = [NPR_eachfile, 1]
       for i_file in 1:parti_file_num
          local p_range = ((i_file - 1) * NPR_eachfile + 1):(i_file * NPR_eachfile)
          local ncfile = NetCDF.open(fileList[i_file], mode=NC_WRITE)

          NetCDF.putvar(ncfile["i"], [parti[n].i for n in p_range], start=start, count=count)
          NetCDF.putvar(ncfile["j"], [parti[n].j for n in p_range], start=start, count=count)
          NetCDF.putvar(ncfile["k"], [parti[n].k for n in p_range], start=start, count=count)
          NetCDF.putvar(ncfile["z"], [parti[n].z for n in p_range], start=start, count=count)
          NetCDF.putvar(ncfile["u"], [parti[n].u for n in p_range], start=start, count=count)
          NetCDF.putvar(ncfile["v"], [parti[n].v for n in p_range], start=start, count=count)
          NetCDF.putvar(ncfile["w"], [parti[n].w for n in p_range], start=start, count=count)
          NetCDF.putvar(ncfile["rho"], [parti[n].rho for n in p_range], start=start, count=count)
          NetCDF.putvar(ncfile["s"], [parti[n].s for n in p_range], start=start, count=count)
          NetCDF.putvar(ncfile["t"], [parti[n].t for n in p_range], start=start, count=count)
          NetCDF.putvar(ncfile["vor"], [parti[n].vor for n in p_range], start=start, count=count)
          NetCDF.putvar(ncfile["shear"], [parti[n].shear for n in p_range], start=start, count=count)
          NetCDF.putvar(ncfile["strain"], [parti[n].strain for n in p_range], start=start, count=count)
          
          NetCDF.sync(ncfile)
         end
      end


  function ini_particles(time)
    #IMPLICIT NONE
    #INTEGER :: i,j,k,time,itmp,npr_eachline
    #REAL(kind=rc_kind) ::  rand,r,theta, x1, y1
    #INTEGER,PARAMETER :: seed = 86456  

    open_parti_files()


    for i in 1:NPR
      parti[i] = Particle(0e0,0e0,0e0,0e0,0e0,0e0,0e0,0e0,0e0,0e0,0e0,0e0,0e0,0e0,0e0,0e0,0e0,0e0,0e0,0e0)
      # set to time?
    end
    
    for i in 1:NPR
      parti[i].i = rc_kind(i*NL/NPR)
      parti[i].j=90e0
      parti[i].k=2e0
    end
 
   end


  function get_parti_vel(time)
    #IMPLICIT NONE
    #INTEGER :: i,j,k,ip,ic,jc,kc,ifc,jfc,kfc,time
    #REAL(kind=rc_kind) ::  dic,djc,dkc,dif,djf,dkf
    #REAL(kind=rc_kind), DIMENSION(    0:NI,0:NJ+1        )        :: uxf
    #REAL(kind=rc_kind), DIMENSION(    0:NI+1,0:NJ        )        :: vyf
    #REAL(kind=rc_kind), DIMENSION(    0:NI+1,0:NJ+1,0:NK )        :: wzf
    #REAL(kind=rc_kind), DIMENSION(      NI,  0:NJ,     NK)        :: vfp
    #REAL(kind=rc_kind), DIMENSION(    0:NI+1,0:NJ,   0:NK+1)      :: vf_ex
    #REAL(kind=rc_kind), DIMENSION(    0:NI,    NJ,     NK)        :: ufp
    #REAL(kind=rc_kind), DIMENSION(    0:NI,  0:NJ+1, 0:NK+1)      :: uf_ex
    #REAL(kind=rc_kind), DIMENSION(      NI,    NJ,   0:NK)        :: wfp
    #REAL(kind=rc_kind), DIMENSION(    0:NI+1,0:NJ+1, 0:NK)        :: wf_ex
   local wzf = OffsetArrays.zeros(rc_kind, 0:NI+1, 0:NJ+1, 0:NK)
   local ufp = OffsetArrays.zeros(rc_kind, 0:NI, 1:NJ, 1:NK)
   local vfp = OffsetArrays.zeros(rc_kind, 1:NI, 0:NJ, 1:NK)
   local wfp = OffsetArrays.zeros(rc_kind, 1:NI, 1:NJ, 0:NK)
   local uf_ex = OffsetArrays.zeros(rc_kind, 0:NI, 0:NJ+1, 0:NK+1)
   local vf_ex = OffsetArrays.zeros(rc_kind, 0:NI+1, 0:NJ, 0:NK+1)
   local wf_ex = OffsetArrays.zeros(rc_kind, 0:NI+1, 0:NJ+1, 0:NK)
    #rearrange the ux and vy to face grids
    #uxf = 0.5d0*(ux(0:NI,:)+ux(1:NI+1,:))
    #vyf = 0.5d0*(vy(:,0:NJ)+vy(:,1:NJ+1))
    @views @. wzf = 0.5e0*(wz[:,:,0:NK] + wz[:,:,1:NK+1])

    #calculate the face velocity
    local k=0
    wfp[:,:,k] = wf[:,:,k]/J2d[1:NI,1:NJ]*wzf[1:NI,1:NJ,k]

    for k in 1:NK
       ufp[:,:,k] = uf[:,:,k]/Jifc[:,:,k]
       vfp[:,:,k] = vf[:,:,k]/Jjfc[:,:,k]
       wfp[:,:,k] = wf[:,:,k]/J2d[1:NI,1:NJ]*wzf[1:NI,1:NJ,k]
    end
    uf_ex .= 0e0
    uf_ex[:,1:NJ,1:NK] = ufp
    # === vertical extrapolation
    uf_ex[:,:,NK+1] = 2*uf_ex[:,:,NK]-uf_ex[:,:,NK-1] # extrapolation

    vf_ex=0d0
    vf_ex[1:NI,:,1:NK] = vfp
    # === zonally periodic
    vf_ex[0,:,:] = vf_ex[NI,:,:]
    vf_ex[NI+1,:,:]=vf_ex[1,:,:]
    # === vertical extrapolation
    vf_ex[:,:,NK+1] = 2*vf_ex[:,:,NK]-vf_ex[:,:,NK-1]

    wf_ex .= 0e0
    wf_ex[1:NI,1:NJ,:] = wfp
    # ===zonal periodic condition
    wf_ex[0,:,:] = wf_ex[NI,:,:]
    wf_ex[NI+1,:,:] = wf_ex[1,:,:]
    # ===
    for ip in 1:NPR
       parti[ip].time= 0e0
       # add time?
       if (parti[ip].j < NJ && parti[ip].j > 0 && parti[ip].k < NK && parti[ip].k > 0)
          #ic, jc, kc, is the integer index of the particle relative to 
          #the grids center. Use these values for variables with the ghost points.
          #ifc, jfc, and kfc is the index relative to the coordinates of grid faces. 
          #Use these values for variables on faces.
          local ic = Int(parti[ip].i+0.5e0)
          local jc = Int(parti[ip].j+0.5e0)
          local kc = Int(parti[ip].k+0.5e0)

          local ifc = Int(parti[ip].i)
          local jfc = Int(parti[ip].j)
          local kfc = Int(parti[ip].k)

          local dif = parti[ip].i - ifc
          local djf = parti[ip].j - jfc
          local dkf = parti[ip].k - kfc

          
          local dic = parti[ip].i - ic + 0.5e0
          local djc = parti[ip].j - jc + 0.5e0
          local dkc = parti[ip].k - kc + 0.5e0
          # calcuate the zonal velocity
          local swap1 = sigma2z(parti[ip].i,parti[ip].j+0.5e0,parti[ip].k+0.5)
          local swap2 = sigma2z(parti[ip].i,parti[ip].j+0.5e0,rc_kind(kc))
          local swap3 = sigma2z(parti[ip].i,parti[ip].j+0.5e0,rc_kind(kc+1))
          local dz = (swap1 - swap2) / ( swap3 - swap2 )

          parti[ip].u = interp_trilinear(dif,djc,dz,uf_ex[ifc:ifc+1,jc:jc+1,kc:kc+1])
          # parti[ip].u = interp_trilinear(dif,djc,dkc,uf_ex[ifc:ifc+1,jc:jc+1,kc:kc+1])

          #calcuate the meridional velocity
          local swap1 = sigma2z(parti[ip].i+0.5e0,parti[ip].j,parti[ip].k+0.5e0)
          local swap2 = sigma2z(parti[ip].i+0.5e0,parti[ip].j,rc_kind(kc))
          local swap3 = sigma2z(parti[ip].i+0.5e0,parti[ip].j,rc_kind(kc+1))
          local dz = (swap1 - swap2) / ( swap3 - swap2 )
          parti[ip].v = interp_trilinear(dic,djf,dz,vf_ex[ic:ic+1,jfc:jfc+1,kc:kc+1])
          #parti[ip].v = interp_trilinear(dic,djf,dkc,vf_ex[ic:ic+1,jfc:jfc+1,kc:kc+1])

          #calcuate the vertical velocity
          local swap1 = sigma2z(parti[ip].i,parti[ip].j,parti[ip].k)
          local swap2 = sigma2z(parti[ip].i,parti[ip].j,rc_kind(kfc))
          local swap3 = sigma2z(parti[ip].i,parti[ip].j,rc_kind(kfc+1))
          local dz = (swap1 - swap2) / ( swap3 - swap2 )
          parti[ip].z = swap1*DL
          parti[ip].w = interp_trilinear(dic,djc,dz,wf_ex[ic:ic+1,jc:jc+1,kfc:kfc+1])
          #parti[ip].w = interp_trilinear(dic,djc,dkf,wf_ex[ic:ic+1,jc:jc+1,kfc:kfc+1])

          #diagnose other properties
          local swap1 = sigma2z(parti[ip].i+0.5e0,parti[ip].j+0.5e0,parti[ip].k+0.5e0)
          local swap2 = sigma2z(parti[ip].i+0.5e0,parti[ip].j+0.5e0,rc_kind(kc))
          local swap3 = sigma2z(parti[ip].i+0.5e0,parti[ip].j+0.5e0,rc_kind(kc+1))
          local dz = (swap1 - swap2) / ( swap3 - swap2 )

          parti(ip).vor = interp_trilinear(dic,djc,dz,vor[ic:ic+1,jc:jc+1,kc:kc+1])
          parti(ip).rho = interp_trilinear(dic,djc,dz,rho[ic:ic+1,jc:jc+1,kc:kc+1])
          parti(ip).shear = interp_trilinear(dic,djc,dz,shear[ic:ic+1,jc:jc+1,kc:kc+1])
          parti(ip).strain = interp_trilinear(dic,djc,dz,strain[ic:ic+1,jc:jc+1,kc:kc+1])
          #interp_trilinear(dic,djc,dkc,vor[ic:ic+1,jc:jc+1,kc:kc+1],parti[ip].vor)
          #interp_trilinear(dic,djc,dkc,s[ic:ic+1,jc:jc+1,kc:kc+1,0],parti[ip].vor)
          #interp_trilinear(dic,djc,dkc,shear[ic:ic+1,jc:jc+1,kc:kc+1],parti[ip].shear)
          #interp_trilinear(dic,djc,dkc,strain[ic:ic+1,jc:jc+1,kc:kc+1],parti[ip].strain)
          #parti[ip].shear = dz
          #parti[ip].strain = dkc
          #parti[ip].r = rho[ic:ic+1,jc:jc+1,kc:kc+1]

       else
          parti[ip].u=0e0
          parti[ip].v=0e0
          parti[ip].w=0e0
       end
      end
    # get the zonal velocity u

   end

  function parti_forward()
    for i in 1:NPR
       parti[i].i = parti[i].i + 0.5e0 * dtf * (3e0 * parti[i].u - parti[i].u0)
       if (parti[i].i > NI)
         parti[i].i = parti[i].i - rc_kind(NI)
       end
       if (parti[i].i <0e0 ) 
         parti[i].i = parti[i].i + rc_kind(NI)
       end

       if (parti[i].j > NJ-1 && parti[i].v > 0)
          parti[i].j = parti[i].j + parti[i].v * dtf / (1e0 + (parti[i].v * dtf)/(rc_kind(NJ) - parti[i].j) )
       elseif (parti[i].j<1 .AND. parti[i].v<0)
          parti[i].j = parti[i].j + parti[i].v * dtf / ( 1e0 - dtf/parti[i].j )
       else
          parti[i].j = parti[i].j + 0.5e0 * dtf * (3e0 * parti[i].v - parti[i].v0)
       end

       if (parti[i].k>NK-1 && parti[i].w>0)
          parti[i].k = parti[i].k + parti[i].w * dtf / (1e0 + (parti[i].w * dtf)/(rc_kind(NK) - parti[i].k) )
       elseif (parti[i].k<1 && parti[i].w<0)
          parti[i].k = parti[i].k + parti[i].w * dtf / ( 1e0 - dtf/parti[i].k ) 
       else
          parti[i].k = parti[i].k + 0.5e0 * dtf * (3e0 * parti[i].w - parti[i].w0)
       end

       #debug part
       if (parti[i].j<0e0 || parti[i].j>NJ || parti[i].k>NK || parti[i].k<0e0 )
          println("particles coordinates are wrong, iPR=",i,"j,k",parti[i].j, " ",parti[i].k)
          exit(1)
       end

       parti[i].u0 = parti[i].u
       parti[i].v0 = parti[i].v
       parti[i].w0 = parti[i].w
      end

   end

  function interp_trilinear(di,dj,dk,var)
    #!== give 8 corner points of a cube, interpolate point values inside of the cube
    #!== di is the distance between the particle to the left face
    #!== dj is the distance between the particle to the southern face
    #!== dk is the distance between the particle and the bottom face
    #IMPLICIT NONE
    #REAL(kind=rc_kind), INTENT(in) :: di,dj,dk
    #REAL(kind=rc_kind), INTENT(in), DIMENSION(    2,    2  ,   2  )      :: var
    #REAL(kind=rc_kind), INTENT(out) :: velp
    #REAL(kind=rc_kind) ::  i1,i2,i3,i4,j1,j2,ti,tj,tk

    # calcuate the Trilinear interpolation
    local i1 = (var[2,1,1] - var[1,1,1])*di + var[1,1,1]
    local i2 = (var[2,1,2] - var[1,1,2])*di + var[1,1,2]
    local i3 = (var[2,2,2] - var[1,2,2])*di + var[1,2,2]
    local i4 = (var[2,2,1] - var[1,2,1])*di + var[1,2,1]

    local j1 = (i3 - i2)*dj + i2
    local j2 = (i4 - i1)*dj + i1

    return (j1 - j2) * dk + j2
  end

  #SUBROUTINE  get_parti_vel_ana()
  #  INTEGER :: ip
  #  DO ip = 1, NPR
  #     parti(ip).u = -1*SIN(pi*parti(ip).i/REAL(NI))*COS(pi*parti(ip).j/REAL(NJ))
  #     parti(ip).v = COS(pi*parti(ip).i/REAL(NI))*SIN(pi*parti(ip).j/REAL(NJ))
  #  ENDDO
  #END SUBROUTINE get_parti_vel_ana


end

end
