function plot_objective(x,e_l,e_u,cost,penalty,epsilon)
% PLOT_OBJECTIVE(x,e_l,e_u,cost,penalty,epsilon)
% Plots the function F(x,e_l,e_r) = max(f(x,e_l),f(x,e_u)) for fixed e_l
% and e_r.
%
% Input:
%       x - 1 by m vector           x(j) corresponds to some fixed hour of
%                                   the j's day
%       e_l - scalar                lower bound on solar radiance
%       e_r - scalar                upper bound on solar radiance
%       cost - scalar               price for an energy unit
%       penalty - function handle   penalty function, should be able to
%                                   act on matrices
%       epsilon - pos. scalar << 1  defines the width of the no-penalty
%                                   interval [x(1-epsilon),x(1+epsilon)]
% Output:
%       no output, function just creates a plot
%% Calculation
x_opt = (penalty(e_u) + cost*e_l)/...                                       % Is an intersection point of two lines:
                    (cost*(1-epsilon) + penalty(1+epsilon));                % penalty*((1-epsilon)x-e_l) and cost*(e_u-(1+epsilon)x)
x1 = e_u/(1+epsilon);                                                       % Is an intersection point of cost*(e_u-(1+epsilon)x) and the X-axis
x2 = e_l/(1-epsilon);                                                       % Is an intersection point of penalty*((1-epsilon)x-e_l) and the X-axis
[y, ~] = medium_objective(x,e_l,e_u,cost,penalty,epsilon);                  % y = F(x,e_l,e_u)
y = -y;
%% Plot
figure, hold on
plot(x,y)                                                                   % plot of F(x,e_l,e_u)
if x1>x2, plot(x_opt,0,'r*'); end                                           % plot the minimizer, if it is a unique one

% no-penalty intervals for some x
x1 = x(y>0); x1 = x1(1:20:end);
y1 = y(y>0); y1 = y1(1:20:end);

plot(x1,y1,'ro')                                                            % central points of no-penalty intervals
plot([x1*(1-epsilon); x1*(1+epsilon)], [y1; y1], 'r--')                     % intervals as lines
plot(x1*(1-epsilon), y1, 'r<')                                              % lower bound
plot(x1*(1+epsilon), y1, 'r>')                                              % upper bound

% confidence interval of p
% if x2<x1
    plot([e_l(1) e_u(1)], [0 0], 'Color', 'black', 'LineWidth', 2)                    % plot the interval [e_l,e_r] as a line
    plot(e_l(1), 0, 'Color', 'black', 'Marker', '<', 'MarkerSize', 10)             % lower bound
    plot(e_u(1), 0, 'Color', 'black', 'Marker', '>', 'MarkerSize', 10)             % upper bound
    plot((e_l(1)+e_u(1))/2, 0, 'ko')                                                  % center
% end

% axis labels
xlabel('x')
ylabel('objective')
hold off