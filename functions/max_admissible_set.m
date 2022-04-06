function [A_set,b_set] = max_admissible_set(A,K,u_lim,x_lim)
    

C = [eye(size(A));
     K];

t=0;

f = [x_lim; u_lim; x_lim; u_lim];
s = length(f);

A_constr = [];
b_constr = [];

out = false;
opts = optimoptions('linprog','Display','off');
while ~out

    A_constr = [A_constr; C*A^t; -C*A^t];
    b_constr = [b_constr; f];
    J = [C*A^(t+1); -C*A^(t+1)];
    
    fprintf('Iteration %i\n',t);
    max_val = zeros(1,s);
    for i=1:s
        [~,fval,exit_opt,~] = linprog(J(i,:),A_constr,b_constr,[],[],[],[],opts);
        fval = -fval;
    
        if exit_opt ==-3
            fval = inf;
        end
        max_val(i) = fval-b_constr(i);
    end
    
    if all(max_val<=0-eps) && exit_opt ~=-3
        out = true;
        A_set = A_constr;
        b_set = b_constr;
    else
        t=t+1;
    end

end
end