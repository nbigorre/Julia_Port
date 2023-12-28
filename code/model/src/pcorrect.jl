function pcorrect(pcorr) 
   @views @inbounds @. p[0:NI+1, 0:NJ+1, 0:NK+1] .+= pcorr[1:NI+2, 1:NJ+2, 1:NK+2]
end                                     
