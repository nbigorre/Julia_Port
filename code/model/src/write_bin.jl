function write_bin(step)
##include "cppdefs.h"
#  USE header
#  IMPLICIT NONE
#
#  INTEGER :: step
#  CHARACTER(LEN=10) :: stepchar 
#  
#
#
#!-------
#! OUTPUT


#if (mod(step,out2d_int) == 0) then
#  save2d(NI+2,NJ+2,rho(:,:,NK),string(dirout,"/op.rho.",lpad(step, 10, "0"),".bin"))
#  save2d(NI+2,NJ+2,h,string(dirout,"/op.h.",lpad(step, 10, "0"),".bin"))
#end
#
#
#if (mod(step,out3d_int) == 0) then
#  save3d(NI+2,NJ+2,NK+2,u(:,:,:,1),string(dirout,"/op.u.",lpad(step, 10, "0"),".bin"))
#  save3d(NI+2,NJ+2,NK+2,v(:,:,:,1),string(dirout,"/op.v.",lpad(step, 10, "0"),".bin"))
#  save3d(NI+2,NJ+2,NK+2,w(:,:,:,1),string(dirout,"/op.w.",lpad(step, 10, "0"),".bin"))
#  save3d(NI+2,NJ+2,NK+2,rho,string(dirout,"/op.rho.",lpad(step, 10, "0"),".bin"))
#  vort(0)
#  save3d(NI+2,NJ+2,NK+2,vor,string(dirout,"/op.vor.",lpad(step, 10, "0"),".bin"))
#end




#-------
# PICKUP

#Reading pickup files, only once, when step is equal to pickup_step
if (pickup_step>-0.5)
  if (step==pickup_step)
    r_pickup(string(dirout,"/op.pickup.",lpad(step, 10, "0"),".bin"))
  end
end


# Writing pickup files at some frequency
if (mod(step,pickup_int) == 0) 
  w_pickup(string(dirout,"/op.pickup.",lpad(step, 10, "0"),".bin"))
end

end
