function [obj, grad] = medium_objective(x,e_l,e_u,cost,penalty,epsilon)
% [obj, grad] = MEDIUM_OBJECTIVE(x,e_l,e_u,cost,penalty,nu)
% Calculates the F(x,e_l,e_u) = max(f(x,e_l),f(x,e_u)) and it's gradient
% (where it does exist). This is a simple model of the 1d cost function for
% the robust optimization approach. 
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
%       obj - n by m matrix         obj(i,j) = max(
%                                 cost*max(0,(1-epsilon)x(i,j)-e_l(j)) +
%                                 penalty(max(0,e_l(j)-(1+epsilon)x(i,j))),
%                                 cost*max(0,(1-epsilon)x(i,j)-e_u(j)) +
%                                 penalty(max(0,e_u(j)-(1+epsilon)x(i,j)))
%                                         )
%       grad - n by m matrix        !!!works only for linear penalties!!!
%                                   pleas see the last part of the code
%% Check the size of x and e
s = size(x);
if size(e_l,2) ~= s(2), error('Sizes of x and e_l do not match'); end
if size(e_u,2) ~= s(2), error('Sizes of x and e_r do not match'); end
%% Calculation
% Objective
obj_caseL = small_objective(x,e_l,cost,penalty,epsilon);
obj_caseU = small_objective(x,e_u,cost,penalty,epsilon);
obj = max(obj_caseL, obj_caseU);

% Gradient
x_opt = ...                                                    % Is an intersection point of two lines:
        (penalty(e_u) + cost*e_l)/(cost*(1-epsilon) + penalty(1+epsilon));  % penalty*((1-epsilon)x-e_l) and cost*(e_u-(1+epsilon)x)
x1 = ones(s(1),1)*e_u/(1+epsilon);                                          % Is an intersection point of cost*(e_u-(1+epsilon)x) and the X-axis
x2 = ones(s(1),1)*e_l/(1-epsilon);                                          % Is an intersection point of penalty*((1-epsilon)x-e_l) and the X-axis
grad = NaN(s);                                                              % Initialization of grad, it will remain NaN at the kinks of function
grad(x<x_opt & x<x1) = cost;                                                % The WC-scenario: we provide more then scheduled
grad(x1<x & x<x2) = 0;                                                      % Interval [e_l,e_r] is covered by the no-penalty interval
grad(x>x_opt & x>x2) = penalty(1);                                          % The WC-scenario: we provide less then scheduled