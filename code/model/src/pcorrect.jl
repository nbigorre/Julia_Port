function pcorrect(pcorr) 
   @views @inbounds @. p[1:NI+2, 1:NJ+2, 1:NK+2] .+= pcorr[1:NI+2, 1:NJ+2, 1:NK+2]
end                                     
