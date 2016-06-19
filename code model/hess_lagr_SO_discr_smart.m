function [ hessian_lagrange ] = hess_lagr_SO_discr_smart( x, E,K,cost,penalty,...
    penalty_grad,penalty_hess,epsilon,P )
%HESS_LAGR_SO_DISCR_SMART outputs the hessian of the lagrange of the
%optimization problem of the smart approach with discretization

% Input:
%       x - 1 by m+3*K*m vector          [x,SOC^1,...,SOC^K,b^{in,1},...,b^{in,K},b^{out,1},...,b^{out,K}] 
%                                        (corresponds to the optimization variable vector in outer opt problem
%                                        start_SO_discr_bat_smart_constr)
%       lambda - structure that contains lambda.ineqnonlin
%       E - L by m vector                E(i,j) is the solar radiance during the
%                                        j^th hour in the i^th sample, L>=K
%       K - scalar                       sample size, i.e. we only use the
%                                        K first rows of E
%       cost - 1 by m vector             price per energy unit (time
%                                        dependent)
%       penalty - function handle        penalty function
%       penalty_grad function handle     derivative of the penalty function
%       epsilon - pos. scalar << 1       width of no-penalty interval
%       P - scalar                       nominal power of PV element
%
% Output:
%       obj - scalar                     weighted objective function value
%       hessian_lagrange - m+3*K*m by m+3*K*m vector  
%                                        hessian of langrange w.r.t. [x,SOC^1,...,SOC^K,b^{in,1},...,b^{in,K},b^{out,1},...,b^{out,K}]
  


[~,~,hessian] = obj_SO_discr_weighted_sum(x,E,K,cost,penalty,penalty_grad,penalty_hess,epsilon,P);
hessian_lagrange=hessian;
end

