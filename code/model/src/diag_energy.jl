function diag_energy(step) 
#---------------------------------------------------                    
# USE header
#                                                                      
# integer :: step 
# REAL(kind=rc_kind) ::  cst,enkin,enkin2
#                                                                  
  local cst= EPS*EPS*delta*delta 

  @views local enkin = sum( Jac[1:NI,1:NJ,1:NK] .* ( u[1:NI,1:NJ,1:NK,0] .^ 2 .+ v[1:NI,1:NJ,1:NK,0] .^ 2 .+ cst .* w[1:NI,1:NJ,1:NK,0] .^ 2  ))
              
  println("#total kinetic energy = ", rpad(step, 10, " "), " ", enkin)
end