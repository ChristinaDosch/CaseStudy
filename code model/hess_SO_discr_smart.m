function hessian = hess_SO_discr_smart(x,x_tilde,penalty_hess,epsilon,P)
% HESSIAN outputs the hessian of the objective w.r.t. [SOC,b^{in]},b^{out}]
% Input: 
%       x - 1 by m vector                schedule, x(i) corresponds to the i^th time step
%       x_tilde - 1 by m vector          x_tilde(i) is the actually provided value at time step i
%       penalty_hess - function handle   hessian of the penalty function
%       epsilon - scalar                 size of the no-penalty interval
%       P - scalar                       nominal power of PV element
% Output:
%       hessian - 3m by 3m matrix        hessian of the objective w.r.t. [SOC,b^{in]},b^{out}]

[yp,gradp,~,~] = smooth_ppart((x - epsilon*P) - x_tilde,1E-3,5);

hessian = zeros(3*T,3*T);

hessian(T+1:2*T,T+1:2*T) = ones(T,1)*(gradp.*penalty_hess(yp).*gradp);
hessian(T+1:2*T,2*T+1:3*T) = ones(T,1)*(- gradp.*penalty_hess(yp).*gradp.*0.95);
hessian(2*T+1:3*T,T+1:2*T) = hessian(T+1:2*T,2*T+1:3*T);
hessian(2*T+1:3*T,2*T+1:3*T) = ones(T,1)*(0.95.*0.95.*gradp.*penalty_hess(yp).*gradp);
end

