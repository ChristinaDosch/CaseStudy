T = 96;
[x_optRO, ~, ~, ~] = start_RO_B_inout(false, 'smooth', 2);
% load('RO plots final/BAT delta = 0.01P, cost = 5, 90% con/output.mat');
x1 = x_optRO(1:T); b_in1 = x_optRO(T+1:2*T); b_out1 = x_optRO(2*T+1:3*T);
% [x2,~,~,~] = start_RO(false, 'smooth');
load('RO plots final/delta = 0.02P, cost = 5, 90% con/output.mat');
x2 = x_optRO(1:T);

[obj1, obj2] = start_comparison(x1,x2,true,'RO with (blue) and without (red) battery',b_in1,b_out1);