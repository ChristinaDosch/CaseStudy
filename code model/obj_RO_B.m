function [obj, grad] = obj_RO_B(xb,e_l,e_u,cost,penalty,penalty_grad,epsilon,P,SmoothOrNonSmooth)
% [obj, grad] = OBJ_RO(xb,e_l,e_u,cost,penalty,epsilon,P,SmoothOrNonSmooth)
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
%       cost - 1 by m vector        cost(j) is a price for an energy unit
%                                   during j's hour
%       penalty - function handle   penalty function, should be able to
%                                   act on matrices
%       epsilon - pos. scalar << 1  defines the width of the no-penalty
%                                   interval [x(1-epsilon),x(1+epsilon)]
%       P - scalar                  nominal power of PV element
%       SmoothOrNonSmooth - logical scalar: true for the smooth version
%                                           false for the nonsmooth version
%
% Output:
%       obj - n by 1 vector         obj(i) = sum_j
%                                               F(x(i,j),e_l(j),e_r(j))
%       grad - n by 1 matrix        !!!works only for linear penalties!!!
%                                   grad(i) = sum_j
%                                               gradF(x(i,j),e_l(j),e_r(j))
%                                   With F from medium_objective
%%
if nargin == 8, SmoothOrNonSmooth = 'nonsmooth'; end
%% Calculation
s = size(xb);
if rem(s(2),2) ~= 0, error('Dimension 2 of xb must be even'), end
x = xb(:,1:round(s(2)/2));
b = xb(:,round(s(2)/2)+1:end);
% whos e_u e_l x b xb s
e_u = e_u - max(b,0);
e_l = e_l + max(-b,0);
% e_u = e_u - b;
% e_l = e_l - b;
[obj, grad] = medium_objective(x,e_l,e_u,cost,penalty,penalty_grad,epsilon,P,SmoothOrNonSmooth);
obj = sum(obj,2);
grad = sum(grad,2);
end

function [obj, grad] = medium_objective(x,e_l,e_u,cost,penalty,penalty_grad,epsilon,P,SmoothOrNonSmooth)
% [obj, grad] = MEDIUM_OBJECTIVE(x,e_l,e_u,cost,penalty,penalty_grad,epsilon,P,SmoothOrNonSmooth)
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
%       cost - 1 by m vector        cost(j) is a price for an energy unit
%                                   during j's hour
%       penalty - function handle   penalty function, should be able to
%                                   act on matrices
%       epsilon - pos. scalar << 1  defines the width of the no-penalty
%                                   interval [x(1-epsilon),x(1+epsilon)]
%       P - scalar                  nominal power of PV element
%       SmoothOrNonSmooth - logical scalar: true for the smooth version
%                                           false for the nonsmooth version
%
% Output:
%       obj - n by m matrix         obj(i,j) = max(
%                                 cost*max(0,(1-epsilon)x(i,j)-e_l(j)) +
%                                 penalty(max(0,e_l(j)-(1+epsilon)x(i,j))),
%                                 cost*max(0,(1-epsilon)x(i,j)-e_u(j)) +
%                                 penalty(max(0,e_u(j)-(1+epsilon)x(i,j)))
%                                         )
%       grad - n by m matrix        !!!works only for linear penalties!!!
%                                   please see the last part of the code

%% Check the size of x and e
s = size(x);
if size(e_l,2) ~= s(2), error('Sizes of x and e_l do not match'); end
if size(e_u,2) ~= s(2), error('Sizes of x and e_r do not match'); end
%% Calculation
% non-smooth
if strcmp(SmoothOrNonSmooth,'nonsmooth'),
    % Objective
    obj_caseL = small_objective(x,e_l,cost,penalty,epsilon,P);
    obj_caseU = small_objective(x,e_u,cost,penalty,epsilon,P);
    obj = max(obj_caseL, obj_caseU);

    % Gradient
    x1 = ones(s(1),1)*(e_l - epsilon*P);
    x2 = ones(s(1),1)*(e_l + epsilon*P);
    cost = ones(s(1),1)*cost;
    grad = NaN(s); % Initialization of grad, it will remain NaN at the kinks of the objective
    I = x<x1;
    grad(I) = cost(I);
    I = (x>x1) & (x<x2);
    grad(I) = 0;
    I = x>x2;
    grad(I) = penalty_grad(x(I));
    
%smooth
elseif strcmp(SmoothOrNonSmooth,'smooth'),
    % Objective
    [obj_caseL, grad_caseL] = small_objective_smooth(x,e_l,cost,penalty,penalty_grad,epsilon,P);
    [obj_caseU, grad_caseU] = small_objective_smooth(x,e_u,cost,penalty,penalty_grad,epsilon,P);
    [obj, gradX, gradY] = smooth_max(obj_caseL, obj_caseU, 0.05);
    
    % Gradient
    grad = gradX.*grad_caseL + gradY.*grad_caseU;
else error('Last argument must be either "smooth" or "nonsmooth"'),
end
end

function obj = small_objective(x,e,cost,penalty,epsilon,P)
% obj = SMALL_OBJECTIVE(x,e,cost,penalty,epsilon,P)
% Calculates a simple 1d cost function f(x,e).
%
% Input:
%       x - n by m matrix           x(i,j) corresponds to the j's hour of
%                                   the i's day
%       e - 1 by m vector           e(j) is the solar radiance during the
%                                   j's hour
%       cost - 1 by m vector        cost(j) is a price for an energy unit
%                                   during j's hour
%       penalty - function handle   penalty function, should be able to act
%                                   on matrices
%       epsilon - pos. scalar << 1  defines the width of the no-penalty
%                                   interval [x(1-epsilon),x(1+epsilon)]
%       P - scalar                  nominal power of PV element
%
% Output:
%       obj - n by m matrix         obj(i,j) =
%                                   cost*max(0,(1-epsilon)x(i,j)-e(j)) +
%                                   penalty(max(0,e(j)-(1+epsilon)x(i,j)))
%% Check the size of x and e
s = size(x);
if size(e,2) ~= s(2), error('Sizes of x and e do not match'); end
%% Calculation
E = ones(s(1),1) * e;                                                       % E is an n by m matrix with e in each column
cost = ones(s(1),1) * cost;
obj = penalty(max(zeros(s), (x - epsilon*P) - E)) + ...
      max(zeros(s), E - (x + epsilon*P)).*cost - E.*cost;                    % just evaluating the formula
end

function [obj, grad] = small_objective_smooth(x,e,cost,penalty,penalty_grad,epsilon,P)
% obj = SMALL_OBJECTIVE_SMOOTH(x,e,cost,penalty,penalty_grad,epsilon,P)
% Calculates a simple 1d cost function f(x,e).
%
% Input:
%       x - n by m matrix           x(i,j) corresponds to the j's hour of
%                                   the i's day
%       e - 1 by m vector           e(j) is the solar radiance during the
%                                   j's hour
%       cost - 1 by m vector        cost(j) is a price for an energy unit
%                                   during j's hour
%       penalty - function handle   penalty function, should be able to act
%                                   on matrices
%       epsilon - pos. scalar << 1  defines the width of the no-penalty
%                                   interval [x(1-epsilon),x(1+epsilon)]
%       P - scalar                  nominal power of PV element
%
% Output:
%       obj - n by m matrix         obj(i,j) =
%                                   cost*max(0,(1-epsilon)x(i,j)-e(j)) +
%                                   penalty(max(0,e(j)-(1+epsilon)x(i,j)))
%% Check the size of x and e
s = size(x);
if size(e,2) ~= s(2), error('Sizes of x and e do not match'); end
%% Calculation
E = ones(s(1),1) * e;                                                       % E is an n by m matrix with e in each column
cost = ones(s(1),1) * cost;
[m1, g1] = smooth_ppart((x - epsilon*P) - E, 0.05);
[m2, g2] = smooth_ppart(E - (x + epsilon*P), 0.05);
obj = penalty(m1) + m2.*cost - E.*cost;                                     % evaluating the formula
grad = penalty_grad(m1).*g1 + g2.*cost;
end