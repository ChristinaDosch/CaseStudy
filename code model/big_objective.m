function [obj, grad] = big_objective(x,e_l,e_u,cost,penalty,epsilon)
% [obj, grad] = BIG_OBJECTIVE(x,p_l,p_u,cost,penalty,nu)
% Calculates the objective function as in (4) from "First approach with RO"
% and its gradient.
% 
% Input:
%       x - n by m matrix           x(i,j) corresponds to the j's hour of
%                                   the i's day
%       e_l - 1 by m vector         e_l(j) is a lower bound on solar
%                                   radiance during the i's hour
%       e_r - 1 by m vector         e_r(j) is an upper bound on solar
%                                   radiance during the i's hour
%       cost - scalar               price for an energy unit
%       penalty - function handle   penalty function, should be able to
%                                   act on matrices
%       epsilon - pos. scalar << 1  defines the width of the no-penalty
%                                   interval [x(1-epsilon),x(1+epsilon)]
%
% Output:
%       obj - n by 1 vector         obj(i) = sum_j
%                                               F(x(i,j),e_l(j),e_r(j))
%       grad - n by 1 matrix        !!!works only for linear penalties!!!
%                                   grad(i) = sum_j
%                                               gradF(x(i,j),e_l(j),e_r(j))
%                                   With F from medium_objective
%% Calculation
[obj, grad] = medium_objective(x,e_l,e_u,cost,penalty,epsilon);
obj = sum(obj,2);
grad = sum(grad,2);