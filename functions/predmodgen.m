function predmod=predmodgen(LTI,dim)

%Prediction matrices generation
%This function computes the prediction matrices to be used in the
%optimization problem

%Prediction matrix from initial state
T=zeros(dim.nx*(dim.N+1),dim.nx);
for k=0:dim.N
    T(k*dim.nx+1:(k+1)*dim.nx,:)=LTI.A^k;
end

%Prediction matrix from input
S=zeros(dim.nx*(dim.N+1),dim.nu*(dim.N));
for k=1:dim.N
    for i=0:k-1
        S(k*dim.nx+1:(k+1)*dim.nx,i*dim.nu+1:(i+1)*dim.nu)=LTI.A^(k-1-i)*LTI.B;
    end
end

% W = zeros(dim.ny,dim.nu*dim.N);
% for i=0:dim.N-1
%     W(1:dim.ny,i*dim.nu+1:(i+1)*dim.nu) = LTI.C*LTI.A^(dim.N-i-1)*LTI.B;
% end

predmod.T=T;
predmod.S=S;
% predmod.W=W;
