using .cppdefs

function evalrho(rhonew, nl)
  @static if (cppdefs.rhoonly)
    evalrho_rho(rhonew, nl)
    #@ccall "./PSOM_LIB.so".evalrho_rho_(pointer(rhonew)::Ptr{rc_kind}, Ref(nl)::Ref{rc_kind})::Cvoid
  else
    evalrho_sT(rhonew, nl)
    #@ccall "./PSOM_LIB.so".evalrho_st_(pointer(rhonew)::Ptr{rc_kind}, Ref(nl)::Ref{Int})::Cvoid
  end
end