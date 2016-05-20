function start_TUMMS(x)
if nargin == 0, x = 0.05 + pdf(makedist('Normal',12,sqrt(12)),linspace(0,24,25)); end
%% Initialize parameters
cost = 1;
penalty = 2; % linear penalty
N = 1 + 24*1; % one hour schedule
n = 1000; % number of the day simulations
t = linspace(0,24,N);

%% Initialize e_l, e_r
pd_e = makedist('Normal',12,sqrt(12));
e = 0.05 + pdf(pd_e,t); % gaussian radiance distribution during the day

mu = 0.01; % variance of the gaussian distribution at each moment
e_l = e - mu;
e_u = e + mu;

pd_E = makedist('Normal',0,mu);
E = random(pd_E,n,N); E = E + ones(n,1)*e;

%% Calculate the cost function
X = ones(n,1)*x;
obj = sum(max(zeros(n,N), X - E)*penalty + ...
          max(zeros(n,N), E - X)*cost, 2);
i_bad = find(obj == max(obj));
i_good = find(obj == min(obj));

%% Plot the solutions and data
%% max/min cost days
figure, hold on
xi = 0:0.1:24;
yi = pchip(t, E(i_good,:), xi); % best weather realization (interpolation)
plot(xi, yi, 'b')               % best weather realization (plot)
yi = pchip(t, E(i_bad ,:), xi); % worst weather realization (interpolation)
plot(xi, yi, 'r')               % worst weather realization (plot)
plot(t, x, 'k*',...             % schedule
     [t; t], [e_l; e_u], 'k+-.')% uncertainty intervals

legend('Tag mit den kleinsten Kosten',...
       'Tag mit den groessten Kosten',...
       'Deine Vorhersage')
xlabel('Zeit'), ylabel('Energie, kWh')
title('Tage mit den kleinsten und den groessten Kosten')
% size
xlim([0 24])
v = axis; % for the positioning of the text boxes
% info-box at the top left corner
text(0.05*v(2),0.98*v(4), ['min Kosten = ', num2str(obj(i_good))], 'Color','blue')
text(0.05*v(2),0.95*v(4), ['max Kosten = ', num2str(obj(i_bad))], 'Color','red')
text(0.05*v(2),0.90*v(4), ['cost = ', num2str(cost)])
text(0.05*v(2),0.87*v(4), ['penalty = ', num2str(penalty)])
hold off
%% all days together
figure, hold on
% plot
plot(t, x, 'r*')                               % schedule
plot(t, E(1:10:end,:), 'ok', 'MarkerSize', 1)  % weather realizations
plot([t; t], [e_l; e_u], 'k+-.')               % uncertainty intervals
% size
xlim([0 24])
v = axis; % for the positioning of the text boxes
% info-box at the top left corner
legend('deine Vorhersage',...
       'simulierte Werte')
xlabel('Zeit'), ylabel('Energie, kWh')
title('Simulierte Werte und die Vorhersage')
text(0.05*v(2),0.98*v(4), ['mean Kosten = ', num2str(mean(obj))])
hold off
