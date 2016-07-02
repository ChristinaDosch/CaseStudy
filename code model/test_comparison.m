T = 96;
[xb_opt1, ~, ~] = start_RO_B_inout(false, 'smooth', 2);
x1 = xb_opt1(1:T); b_in1 = xb_opt1(T+1:2*T); b_out1 = xb_opt1(2*T+1:3*T);
[x2,~,~] = start_RO(false, 'smooth');
[obj1, obj2] = start_RO_comparison(x1,x2,true,'RO with (red) and without (blue) battery, tight ramping',b_in1,b_out1);