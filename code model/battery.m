function [x_tilde] = battery(e, x, C)
%BATTERY computes x_tilde for a given realization e and a schedule x
% 
% Input:
%       x - 1 by T vector           x(i) corresponds to the i^th hour
%       e - 1 by T vector           e(i) is the solar radiance during the
%                                   i^th hour
%       C - scalar                  capacity of the battery
%
% Output:
%       x_tilde - 1 by T vector     what is really fed in the grid, based
%                                   on the schedule and the given realization

%% Check the size of x and e
T = size(x,2);
if size(e,2) ~= s, error('Sizes of x and e do not match'); end

%% Calculation
% In the following, I drop the "tilde" in the variable names for reasons of
% clarity

b_in = zeros(1,T); % \tilde{b_in}
b_out = zeros(1,T); % \tilde{b_out}
x_tilde = zeros(1,T); % \tilde{x}, to be returned
SOC = zeros(1,T+1); % state of charge ((T+1)-dim since we need a dummy value for easier computation using the for-loop)

SOC_0 = 0.25; % 25 percent state of charge(SOC) at day break
SOC(1) = SOC_0;

for i = 1:T                    
    b_in(i) = min(max(0,e(i)-x(i)),(0.95-SOC(i))*C);               % (11) in the current documentation
    b_out(i) = min(max(0,1/0.95*(x(i)-e(i))),(SOC(i)-0.1)*C);      % (12)
    x_tilde(i) = e(i) - b_in(i) + 0.95*b_out(i);                   % (13)
    SOC(i+1) = SOC(i) + 1/C*(0.95*b_in(i) - b_out(i));             % (10) 
end

% TO DO: Find a way to express this for loop by vector and matrix calculations

end

