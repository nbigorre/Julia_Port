using .header
using .fortVar

function cfdiv(maxdiv::Ref{rc_kind})
    maxdiv[] = 0e0
    local cxf = @fortGetArr("cxf", rc_kind, (NI+1, NJ, NK))
    local cyf = @fortGetArr("cyf", rc_kind, (NI, NJ+1, NK))
    local czf = @fortGetArr("czf", rc_kind, (NI, NJ, NK+1))
    for k in 1:NK
        for j in 1:NJ
            for i in 1:NI
                local cxfdx = (cxf[i+1,j,k] - cxf[i,j,k])
                local cyfdy = (cyf[i,j+1,k] - cyf[i,j,k])
                local czfdz = (czf[i,j,k+1] - czf[i,j,k])
                                        
                local div= abs(cxfdx+ cyfdy + czfdz) 
                if (div > maxdiv[])
                    maxdiv[] = div
                end
            end
        end
    end
end