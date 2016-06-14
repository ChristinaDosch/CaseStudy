function [SOC,b_in,b_out] = battery_smart(e,x_plan,T,P,cost,penalty,penalty_grad,epsilon,...
    A_smart, b_smart, SOC_min, SOC_max)
%BATTERY_SMART computes optimal SOC, b_in, b_out for a given realization e and a schedule x
% 
% Input:
%       x_plan - 1 by T vector      x(i) corresponds to the i^th hour
%       e - 1 by T vector           e(i) is the solar radiance during the
%                                   i^th hour
%       C - scalar                  capacity of the battery
%       SOC_0 - scalar              state of charge(SOC) at day break
%                                   (usually SOC_0 = 0.25)
%
% Output:
%       SOC - 1 by T vector     
%       b_in - 1 by T vector
%       b_out - 1 by T vector    
%       
%  
%%
%[x_min, x_max, delta,~,~,~,~, A_smart, b_smart, SOC_min, SOC_max] = init_constraints(T,P,C,SOC_0);
x_tilde = @(x) e + 0.95*x(2*T+1:3*T) - x(T+1:2*T);
F = @(x) obj_SO_discr(x_plan,x_tilde(x),cost,penalty,penalty_grad,epsilon,P);

x0 = zeros(1,3*T);
options=optimoptions('fmincon', 'MaxFunEvals', 30000); %'SpecifyObjectiveGradient',true);
[x_opt] = fmincon(F, x0, A_smart(2*T-1:4*T-2,T+1:4*T), b_smart(2*T-1:4*T-2),[],[],...
    [SOC_min*ones(1,T),0*ones(1,(2*T))],[SOC_max*ones(1,T),2*P*ones(1,(2*T))], [], options);

SOC = x_opt(1:T);
b_in = x_opt(1*T+1:2*T);
b_out = x_opt(2*T+1:3*T);
end

