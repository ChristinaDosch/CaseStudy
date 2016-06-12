function [obj,grad] = obj_SO_discr(x,e,cost,penalty,epsilon,P)
% obj = REVENUE(x,e,cost,penalty,P)
% Calculates the objective function -F(x,E) as in 1.5 Objective function
% Since F(x,E) = \sum_{i=1}^T F^{(i)}(x_i,E_i) (see also 1.5), obj_SO_discr
% calculates the revenue values F^{(i)}(x_i,E_i) for all time steps 
% and returns the sum of all values
% 
% Input:
%       x - n by m matrix           x(i,j) corresponds to the j's hour of
%                                   the i's day
%       e - 1 by m vector           e(j) is the solar radiance during the
%                                   i's hour
%       cost - scalar               price for an energy unit
%       penalty - function handle   penalty function, should be able to
%                                   act on matrices
%       epsilon - pos. scalar << 1  defines the width of the no-penalty
%                                   interval [x(1-epsilon),x(1+epsilon)]
%       P - scalar                  nominal power of PV element
%
% Output:
%       obj - n by 1 vector         obj(i) = sum_j
%                                               F(x(i,j),e(j))

%% Check the size of x and e
s = size(x);
if size(e,2) ~= s(2), error('Sizes of x and e do not match'); end
%% Calculation
E = ones(s(1),1) * e;                                                       % E is an n by m matrix with e in each column
[yp,gradp,pp]=smooth_ppart((x - epsilon*P) - E,1E-3,5);
[yc,gradc,pc]=smooth_ppart(E - (x + epsilon*P),1E-3,5);
obj = penalty(yp) + yc*cost - E*cost;  

grad = 2*yp.*gradp + gradc*cost ;

%obj = penalty(max(zeros(s), (x - epsilon*P) - E)) + ...
%      max(zeros(s), E - (x + epsilon*P))*cost - E*cost;                     % just evaluating the formula

obj = sum(obj,2); % sum of all single revenue values