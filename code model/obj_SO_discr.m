function [obj,grad_SOC_b,grad_x,hessian] = obj_SO_discr(x,e,cost,penalty,penalty_grad,epsilon,P,penalty_hess)
%function [obj,grad] = obj_SO_discr(x,e,cost,penalty,penalty_grad,P)
% Calculates the objective function -F(x,E) as in 1.5 Objective function
% Since F(x,E) = \sum_{i=1}^T F^{(i)}(x_i,E_i) (see also 1.5), obj_SO_discr
% calculates the revenue values F^{(i)}(x_i,E_i) for all time steps 
% and returns the sum of all values
% 
% Input:
%       x - n by m matrix                x(i,j) corresponds to the j's hour of
%                                        the i's day
%       e - 1 by m vector                e(j) is the solar radiance during the
%                                        i's hour
%       cost - 1 by m vector             price for an energy unit
%       penalty - function handle        penalty function, should be able to
%                                        act on matrices
%       penalty_grad function handle     derivative of the penalty function
%       epsilon - pos. scalar << 1       defines the width of the no-penalty
%                                        interval [x(1-epsilon),x(1+epsilon)]
%       P - scalar                       nominal power of PV element
%       var - text                       variables w.r.t. which the
%                                        gradient is calculated
% Output:
%       obj - n by 1 vector              obj(i) = sum_j
%                                               F(x(i,j),e(j))
%       grad_SOC_b - n by 3m matrix      grad w.r.t. [SOC,b^{in},b^{out}]
%       grad_x - n by m matrix           grad w.r.t. x
%       hessian - 4m by 4m matrix        hessian of the objective w.r.t. [x,SOC,b^{in]},b^{out}]
%% Check the size of x and e
s = size(x); T = s(2); 
if size(e,2) ~= s(2), error('Sizes of x and e do not match'); end
if s(1) ~= 1, error('x schedules for than one day'); end
%% Computation of the objective
E = ones(s(1),1) * e;                        % E is an n by m matrix with e in each row
Cost = ones(s(1),1) * cost;                  % C is an n by m matrix with cost in each row
[yp,gradp,hessp,~] = smooth_ppart((x - epsilon*P) - E,1E-3,5);
[yc,gradc,hessc,~] = smooth_ppart(E - (x + epsilon*P),1E-3,5);

obj = penalty(yp) + yc.*Cost - E.*Cost;  
obj = sum(obj,2); % sum of all single revenue values

%% Computation of gradient w.r.t x
grad_x = penalty_grad(yp).*gradp - gradc.*Cost; 

%% Computation of gradient w.r.t [SOC,b^in,b^out]
Cost = ones(s(1),1) * [cost, cost, cost]; % C is an n by 3m matrix with three times cost in each row
grad_x_tilde = [zeros(s(1),T), -1 * ones(s(1),T), 0.95 * ones(s(1),T)];
grad_SOC_b = - penalty_grad([yp,yp,yp]).*[gradp, gradp, gradp].*grad_x_tilde + grad_x_tilde.*Cost.*[gradc, gradc, gradc] - grad_x_tilde.*Cost ;

%% Computation of the hessian w.r.t. [x,SOC,b^in,b^out]
hessian = zeros(4*T,4*T);
% hessian(1:T,1:T) = ones(T,1)*(penalty_hess(yp).*gradp.*gradp + penalty_grad(yp).*hessp + cost.*hessp);
% hessian(1:T,2*T+1:3*T) = ones(T,1)*(penalty_hess(yp).*gradp.*1.*gradp + penalty_grad(yp).*hessp + cost.*hessp);
% hessian(1:T,3*T+1:4*T) = ones(T,1)*(penalty_hess(yp).*gradp.*(-0.95).*gradp + penalty_grad(yp).*hessp + cost.*hessp);
% hessian(2*T+1:3*T,1:T) = hessian(1:T,2*T+1:3*T);
% hessian(3*T+1:4*T,1:T) = hessian(1:T,3*T+1:4*T);
% hessian(2*T+1:3*T,2*T+1:3*T) = ones(T,1)*(gradp.*penalty_hess(yp).*gradp);
% hessian(2*T+1:3*T,3*T+1:4*T) = ones(T,1)*(- gradp.*penalty_hess(yp).*gradp.*0.95);
% hessian(3*T+1:4*T,2*T+1:3*T) = hessian(T+1:2*T,2*T+1:3*T);
% hessian(3*T+1:4*T,3*T+1:4*T) = ones(T,1)*(0.95.*0.95.*gradp.*penalty_hess(yp).*gradp);
%% Alter Code:
%obj = penalty(max(zeros(s), (x - epsilon*P) - E)) + ...
%   max(zeros(s), E - (x + epsilon*P))*cost - E*cost;                     % just evaluating the formula
end

