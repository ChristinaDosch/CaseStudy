function [hessian] = hess_SO_discr_smart(x_plan,x_tilde,penalty_hess,epsilon,P)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[yp,gradp,~] = smooth_ppart((x_plan - epsilon*P) - x_tilde,1E-3,5);

hessian = zeros(3*T,3*T);

hessian(T+1:2*T,T+1:2*T) = gradp.*penalty_hess(yp).*gradp;
hessian(T+1:2*T,2*T+1:3*T) = - gradp.*penalty_hess(yp).*gradp.*0.95;
hessian(2*T+1:3*T,T+1:2*T) = hessian(T+1:2*T,2*T+1:3*T);
hessian(2*T+1:3*T,2*T+1:3*T) = 0.95.*0.95.*gradp.*penalty_hess(yp).*gradp;
end

