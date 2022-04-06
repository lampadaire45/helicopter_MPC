function [H,h]=costgen_tracking(predmod,weight,dim)
if dim.Nref == 0
    dim.Nref = 1;
end

Qbar=blkdiag(kron(eye(dim.N),weight.Q),weight.beta*weight.P);
Rbar=kron(eye(dim.N),weight.R);
H=predmod.S'*Qbar*predmod.S+Rbar;   
hx0=predmod.S'*Qbar*predmod.T;
hxref=-predmod.S'*Qbar*kron(ones(dim.N+1,dim.Nref),eye(dim.nx));
huref=-Rbar*kron(ones(dim.N,dim.Nref),eye(dim.nu));
h=[hx0 hxref huref];
 
end