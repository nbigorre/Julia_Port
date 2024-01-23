using .fortVar

function cfdiv(maxdiv::Ref{rc_kind})
    maxdiv[] = 0e0
    for k in 1:NK
        for j in 1:NJ
            for i in 1:NI
                local cxfdx = (cxf[i,j,k] - cxf[i-1,j,k])
                local cyfdy = (cyf[i,j,k] - cyf[i,j-1,k])
                local czfdz = (czf[i,j,k] - czf[i,j,k-1])
                                        
                local div= abs(cxfdx+ cyfdy + czfdz) 
                if (div > maxdiv[])
                    maxdiv[] = div
                end
            end
        end
    end
end