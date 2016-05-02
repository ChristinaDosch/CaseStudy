function [obj,grad] = obj_SO_closed_form(x,H1,H2,cost,penalty,epsilon)
%calculates the objective value using the closed-form formula for the
%expected value E[F(,x,E)] as in (3) in section "1.2 Ansatz 1: Stochastic
%Optimization" in objective.tex 
%(however, we consider a min problem here, i.e. the formula is multiplied by -1)

int1 = integral(H1,0,x);
int2 = integral(H2,0,x);

obj = - cost*(1+epsilon)*x + cost*(1+epsilon)*int1 + penalty((1-epsilon)*int2);

grad = -cost*(1+epsilon) + cost*(1+epsilon)*H1(x) + penalty((1-epsilon)*H2(x));

end

