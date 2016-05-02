% 1.) Initialization of parameters, objective function and constraints.
% 2.) Performing robust optimization.
% 3.) Plotting solutions and data
%% Initialize parameters
cost = 1;
epsilon = 0.05;
penalty = @(x) x*2; % linear penalty
N = 1 + 24*1; % one hour schedule
t = linspace(0,24,N);

%% Initialize e_l, e_r for RO
pd = makedist('Normal');
e = 0.2 + pdf(pd,(t-12)/sqrt(12)); % gaussian radiance distribution during the day
mu = 0.1; % uncertanty interval width parameter
e_l = e*(1-mu);
e_u = e*(1+mu);

%% Initialize optimization model for RO
x0 = e; % Starting guess for pattern search
objfct = @(x) big_objective(x,e_l,e_u,cost,penalty,epsilon); % Objective function for RO

% constraints, see "First approach with RO" for explanation
x_min = 0;  % tba
x_max = 0.55; % tba
delta = 0.04; % tba
B = [-eye(N-1) zeros(N-1,1)] + [zeros(N-1,1) eye(N-1)];
A = [B; -B];
b = ones(2*(N-1),1)*delta;

%% Performing optimization
% [x_opt, obj_opt] = ga(objfct,N,A,b,[],[],0,x_max); % genetic algorithm
tic
[x_opt, obj_opt] = patternsearch(objfct,x0,A,b,[],[],x_min*ones(1,N),x_max*ones(1,N)); % pattern search
toc

%% Plot the solutions and data
figure, hold on
plot(t,x_opt,'*r',... % solution computed by ga or patternsearch
     t,(1+epsilon)*x_opt,'^r',...
     t,(1-epsilon)*x_opt,'vr',...
     [t(1) t(end)], [x_max x_max], 'k--',... % x_max
     [t; t], [e_l; e_u], 'k+-.',... % uncertainty intervals
     [t(1) t(end)], [x_min x_min], 'k--',... % x_min
     t,e,'ko') % centers of uncertainty intervals
legend('calculated opt. sol.',...
       'upper no-penalty bound',...
       'lower no-penalty bound',...
       'x_{max}, x_{min}',...
       'uncertainty intervals')
xlabel('time'), ylabel('energy, kWh')
% size
xlim([0 24])
v = axis;
ylim([0 v(4)]);
% info-box at the top left corner
text(0.05*v(2),0.98*v(4),['optimal value = ', num2str(obj_opt)])
text(0.05*v(2),0.93*v(4),['cost = ', num2str(cost)])
text(0.05*v(2),0.90*v(4),['penalty = ', func2str(penalty)])
text(0.05*v(2),0.87*v(4),['[x_{min} x_{max}] = ', '[', num2str(x_min), ' ', num2str(x_max), ']'])
text(0.05*v(2),0.84*v(4),['\Delta = ', num2str(delta)])
hold off
