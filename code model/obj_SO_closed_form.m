function [obj,grad] = obj_SO_closed_form(x,H,exp,cost,penalty,epsilon,P)
%calculates the objective value using the closed-form formula for the
%expected value E[F(x,E)] as in (14) in section "3.1 Ansatz 1: Stochastic
%Optimization"
%(however, we consider a min problem here, i.e. the formula is multiplied by -1)
%
%It can be used for a general number of time steps (i.e. T>1 is possible)
%and various distributions for the E_i

% Input:
%       x: T by 1 vector                 x(i) corresponds to the i-th hour of the day
%       H: cell of function handles      distribution functions of E_i
%       exp: T by 1 vector               vector with expected values for E       
%       cost: scalar                     price for an energy unit
%       penalty: function handle         penalty function
%       epsilon: pos. scalar << 1        defines the width of the no-penalty
%                                        interval [x-epsilon*P,x+epsilon*P]
%       P: scalar                        nominal power of PV element
%
% Output:
%       obj: scalar                 obj = sum_i F^{(i)}(x(i),e(i))
%       grad: scalar

T = size(x,1);
int1 = zeros(T,1); % vector containing the intervals int_0^x(i) H_i(z+epsilon*P)dz
int2 = zeros(T,1);

H_help=zeros(T,1);  % vector containing H_i(epsilon*P)
H1_help=zeros(T,1); % vector containing H_i(x(i)+epsilon*P)
H2_help=zeros(T,1); % vector containing H_i(x(i)-epsilon*P)

for i=1:T
H1=@(z) H{i}(z+epsilon*P);
H2=@(z) H{i}(z-epsilon*P);  

int1(i) = integral(H1,0,x(i));
int2(i) = integral(H2,0,x(i));

H_help(i)=H1(0);     % = H{i}(epsilon*P)
H1_help(i)=H1(x(i)); % = H{i}(x(i)+epsilon*P)
H2_help(i)=H2(x(i)); % = H{i}(x(i)-epsilon*P)
    
end

% calculate objective as given by closed-form formula:
obj = - cost*epsilon*P*ones(T,1) - cost*(H_help.*(exp - epsilon*P*ones(T,1))) - cost*x + cost*int1 + penalty(int2);
% calculate gradient w.r.t x:
grad = -cost*ones(T,1) + cost*H1_help + penalty(H2_help);

end

