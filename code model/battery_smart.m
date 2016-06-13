function [SOC,b_in,b_out] = battery_smart(e,x_plan,C,SOC_0,T,P,cost,penalty,penalty_grad,epsilon)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[x_min, x_max, delta,~,~,~,~, A_smart, b_smart, SOC_min, SOC_max] = init_constraints(T,P,C,SOC_0);
x_tilde = @(x) e + 0.95*x(3*T+1:4*T) - x(2*T+1:3*T);
F = @(x) obj_SO_discr(x_plan
,x_tilde(x),cost,penalty,penalty_grad,epsilon,P);

x0 = zeros(1,4*T);
[x_opt, obj_opt] = fmincon(F, x0, A_smart, b_smart,[],[],...
    [x_min*ones(1,T), SOC_min*ones(1,T),0*ones(1,(2*T))],[x_max*ones(1,T), SOC_max*ones(1,T),2*P*ones(1,(2*T))]);

SOC = x_opt(T+1:2*T);
b_in = x_opt(2*T+1:3*T);
b_out = x_opt(3*T+1:4*T);
end

