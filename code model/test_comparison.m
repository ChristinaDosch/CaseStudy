T = 96;
% [x_optRO, ~, ~, ~] = start_RO_B_inout(false, 'smooth', 2);
load('RO plots final/delta = 0.02P, cost = 5, 90% con/output.mat');
x2 = x_optRO;

% x1 = x_optRO(1:T); b_in1 = x_optRO(T+1:2*T); b_out1 = x_optRO(2*T+1:3*T);
% [x2,~,~,~] = start_RO(false, 'smooth');
% load('SO plots final/SO_x_opt_for_comparison_no_battery.mat');
[x_optRO, ~, ~, ~] = start_RO_B_inout(false, 'smooth', 2);
x1 = x_optRO(1:T); b_in1 = x_optRO(T+1:2*T); b_out1 = x_optRO(2*T+1:3*T);

[obj1, obj2] = start_comparison(x1,x2,true,'RO realistic pricing', b_in1, b_out1);