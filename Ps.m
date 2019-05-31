function [ Ps] = Ps(taufall,PsPeak, time )

Ps=PsPeak.*exp((-1/taufall).*(time)); 

end

