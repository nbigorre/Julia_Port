

function evalrho_sT(rhonew, nl2)

   for k in 0:NK+1
      for j in 0:NJ+1
         for i in 0:NI+1
            local se = s[i+1, j+1, k+1, nl2+1]
            local Te = T[i+1, j+1, k+1, nl2+1]
            rhonew[i+1, j+1, k+1] = potdens(se, Te)
         end
      end
   end

end