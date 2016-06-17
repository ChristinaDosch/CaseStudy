function hessian = hessianfcn(x,lambda)
% HESSIAN returns hessian of the Lagrangean function
% hessian = hessian of objective + sum lambda_i * hessian(inequality
% constraints) + sum lambda_i * hessian(equality constraints)
% Input: 
%        x - 1 by T+3*KT vector         =[x,SOC^1,...,SOC^K,b^{in,1},...,b^{in,K},b^{out,1},...,b^{out,K}]
%        lambda - 1 by T+3*KT vector    Lagrange multipliers
%
%
%

s = size(x,2);
%[obj,grad_SOC_b,grad_x,hessian] = obj_SO_discr(x(1:T),e,cost,penalty,penalty_grad,epsilon,P,penalty_hess)
hessian = zeros(s,s);

% TO DO:
% As soon as computation of hessian of the objective is done, we can
% continue with the hessians of the constraints in order to obtain the
% hessian of the Lagrangean function
end

