function [y, grad, hess, p] = smooth_ppart(x,epsilon,N)
% SMOOTH_PART(x,epsilon) berechnet den positiven Teil von x, falls
% |x|>=epsilon, oder epsilon/4 + 0.5*x + 0.25/epsilon*x^2, falls
% |x|<epsilon. Das entspricht der Gl?ttung der Funktion max(x,0) in dem
% Intervall [-epsilon epsilon]. Here a quadratic polynomial interpolation
% between the points -epsilon and epsilon is used.
% 
% Input:
%       x - n by m matrix           
%       epsilon - pos. scalar << 1            
%       N - natural number > 5      can be used to produce twice
%                                   differentiable approximation
%
% Output:
%       y - n by m matrix           value of the smooth function at x
%       grad - n by m matrix        value of the derivative at x
%       hess - n by m matrix        value of the 2nd derivative at x 
%   !!!! hess is only correctly implemented for N=2 !!!!
%       p - 1 by 3 or
%           1 by N+1 vector         vector of coefficients from the
%                                   interpolating polynomial

s = size(x);

%% the first case is only relevant, if a twice differentiable function is needed
if nargin == 3,
    A = [epsilon.^(0:N);...
        (-epsilon).^(0:N);...
        (0:N).*[1 epsilon.^(0:N-1)];...
        (0:N).*[1 (-epsilon).^(0:N-1)];...
        (-1:N-1).*(0:N).*[1 1 epsilon.^(0:N-2)];...
        (-1:N-1).*(0:N).*[1 1 (-epsilon).^(0:N-2)]
        ];
    A = A(:,1:end-1);
    b = [epsilon; 0; 1; 0; 0; 0] - epsilon^N;
    p = linsolve(A,b);
    p = [1; p(end:-1:1)]';

elseif nargin == 2, N = 2; p = [0.25/epsilon 0.5 epsilon/4]; % quadratic interpolation with precomputed polynomial
end

y = zeros(s);
grad = zeros(s);
hess = zeros(s);

%% |x|<eps
I = (x > -epsilon) & (x < epsilon);
y(I) = polyval(p,x(I));
grad(I) = polyval((N:-1:1).*p(1:end-1), x(I));
hess(I) = polyval(N.*p(1:end-2), x(I));
%% x<-eps
I = (x <= -epsilon);
y(I) = 0;
grad(I) = 0;
hess(I) = 0;
%% x>eps
I = (x >= epsilon);
y(I) = x(I);
grad(I) = 1;
hess(I) = 0;
