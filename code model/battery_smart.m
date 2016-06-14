function x_tilde_opt = battery_smart(e,x_plan,T,P,cost,penalty,penalty_grad,epsilon,...
    A_smart, b_smart, SOC_min, SOC_max)
%BATTERY_SMART computes optimal SOC, b_in, b_out for a given realization e and a given schedule x
% 
% Input:
%       e - 1 by T vector                e(i) is the solar radiance during the i^th hour
%       x_plan - 1 by T vector           given schedule, x(i) corresponds to the i^th hour
%       T - scalar                       number of time steps
%       P - scalar                       nominal power
%       cost - 1 by T                    cost(j) is the price for an energy unit during the j^th hour
%       penalty - function handle        penalty function
%       penalty_grad - function handle   derivative of the penalty function
%       epsilon                          size of no-penalty interval
%       A_smart, b_smart                 contain constraints with battery
%       SOC_min, SOC_max                 lower and upper bounds for SOC
%
% Output:
%       SOC - 1 by T vector     
%       b_in - 1 by T vector
%       b_out - 1 by T vector    
%       
%  
%%
if size(e,2) ~= size(x_plan,2), error('Sizes of x and e do not match'); end
if size(e,2) ~= T, error('Size of x and e is not T'); end

% x \in R^3T is going to be the optimization variable in this function. x = [SOC,b^in,b^out]

x_tilde = @(x) e + 0.95*x(2*T+1:3*T) - x(T+1:2*T);
F = @(x) obj_SO_discr(x_plan,x_tilde(x),cost,penalty,penalty_grad,epsilon,P); % for a given e and x_plan, the objective is a function of b^in and b^out

x0 = zeros(1,3*T);
options = optimoptions('fmincon', 'MaxFunEvals', 30000); %'SpecifyObjectiveGradient',true);
[x_opt] = fmincon(F, x0, A_smart(2*T-1:4*T-2,T+1:4*T), b_smart(2*T-1:4*T-2),[],[],...
    [SOC_min*ones(1,T),0*ones(1,(2*T))],[SOC_max*ones(1,T),2*P*ones(1,(2*T))], [], options); % calculates optimal SOC, b^in, b^out

SOC = x_opt(1:T);
b_in = x_opt(1*T+1:2*T);
b_out = x_opt(2*T+1:3*T);

x_tilde_opt = e + 0.95*b_out - b_in;
end

