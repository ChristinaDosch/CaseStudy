function [obj,grad] = obj_SO_closed_form(x,H,exp,cost,penalty,epsilon,P)
%calculates the objective value using the closed-form formula for the
%expected value E[F(x,E)] as in (14) in section "3.1 Ansatz 1: Stochastic
%Optimization"
%(however, we consider a min problem here, i.e. the formula is multiplied by -1)

% Input:
%       x: T by 1 vector            x(i) corresponds to the i-th hour of the day
%       H: function handle          distribution function of E (should become multivariate)
%       exp: T by 1 vector          vector with expected values for E       
%       cost: scalar                price for an energy unit
%       penalty: function handle    penalty function
%       epsilon: pos. scalar << 1   defines the width of the no-penalty
%                                   interval [x-epsilon*P,x+epsilon*P]
%       P: scalar                   nominal power of PV element
%
% Output:
%       obj: scalar                 obj = sum_i F^{(i)}(x(i),e(i))
%       grad: scalar

T = size(exp,1);
int1 = zeros(T,1);
int2 = zeros(T,1);

H1=@(z) H(z+epsilon*P);
H2=@(z) H(z-epsilon*P);
for i=1:T
int1(i) = integral(H1,0,x(i));
int2(i) = integral(H2,0,x(i));
end

obj = - cost*epsilon*P*ones(T,1) - cost*H(epsilon*P)*(exp - epsilon*P*ones(T,1)) - cost*x + cost*int1 + penalty(int2);

grad = -cost*ones(T,1) + cost*H1(x)*ones(T,1) + penalty(H2(x));

end

