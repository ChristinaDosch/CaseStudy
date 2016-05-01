function obj = small_objective(x,e,cost,penalty,epsilon)
% obj = SMALL_OBJECTIVE(x,e,cost,penalty,epsilon)
% Calculates a simple 1d cost function f(x,e).
%
% Input:
%       x - n by m matrix           x(i,j) corresponds to the j's hour of
%                                   the i's day
%       e - 1 by m vector           e(j) is the solar radiance during the
%                                   i's hour
%       cost - scalar               price for an energy unit
%       penalty - function handle   penalty function, should be able to act
%                                   on matrices
%       epsilon - pos. scalar << 1  defines the width of the no-penalty
%                                   interval [x(1-epsilon),x(1+epsilon)]
%
% Output:
%       obj - n by m matrix         obj(i,j) =
%                                   cost*max(0,(1-epsilon)x(i,j)-e(j)) +
%                                   penalty(max(0,e(j)-(1+epsilon)x(i,j)))
%% Check the size of x and e
s = size(x);
if size(e,2) ~= s(2), error('Sizes of x and p do not match'); end
%% Calculation
E = ones(s(1),1) * e;                                                       % E is an n by m matrix with e in each column
obj = max(zeros(s), x*(1-epsilon) - E) * cost + ...
      penalty(max(zeros(s), E - x*(1+epsilon)));                            % just evaluating the formula