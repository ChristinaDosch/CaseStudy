function obj = revenue(x,e,cost,penalty,epsilon,P)
% obj = REVENUE(x,e,cost,penalty,P)
% Calculates the objective function F(x,E) as in 1.5 Objective function
% Since F(x,E) = \sum_{i=1}^T F^{(i)}(x_i,E_i) (see also 1.5), revenue
% calls small_objective and returns the sum of all values
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
%% Calculation
obj = small_objective(x,e,cost,penalty,epsilon,P);
obj = sum(obj,2);