

function evalrho_rho(rhonew, n)

   for k in 0:NK+1
      for j in 0:NJ+1
         for i in 0:NI+1
            rhonew[i+1, j+1, k+1] = s[i+1, j+1, k+1, n+1]
         end
      end
   end

end