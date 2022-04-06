function[xr,ur]=optimalss(LTI,dim,dtilde)

eqconstraints.A=[eye(dim.nx)-LTI.A -LTI.B; LTI.C zeros(dim.ny,dim.nu)];
eqconstraints.b=[zeros(dim.nx,1); LTI.yref-LTI.Cd*dtilde];

H=eye(dim.nu+dim.nx);
h=zeros(dim.nx+dim.nu,1);


options1 = optimoptions(@quadprog); 
% options1.Optimalityoptions1 = optimoptions(@quadprog); Tolerance=1e-20;
options1.ConstraintTolerance=1.0000e-15;
options1.Display='off';

xur=quadprog(H,h,[],[],eqconstraints.A,eqconstraints.b,[],[],[],options1);
xr=xur(1:dim.nx);
ur=xur(dim.nx+1:end);

end