function [obj,grad] = obj_SO_closed_form(x,H,exp,cost,penalty,epsilon,P)
%calculates the objective value using the closed-form formula for the
%expected value E[F(x,E)] as in (14) in section "3.1 Ansatz 1: Stochastic
%Optimization"
%(however, we consider a min problem here, i.e. the formula is multiplied by -1)

H1=@(z) H(z+epsilon*P);
H2=@(z) H(z-epsilon*P);
int1 = integral(H1,0,x);
int2 = integral(H2,0,x);

obj = - cost*epsilon*P - cost*H(epsilon*P)*(exp - epsilon*P) - cost*x + cost*int1 + penalty(int2);

grad = -cost + cost*H1(x) + penalty(H2(x));

end

