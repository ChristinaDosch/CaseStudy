function start_TUMMS
%% Initialize parameters
P = 20;
N = 1 + 24*1; % one hour schedule
t = linspace(0,24,N);
x0 = 1;
str = load('Optimal x');
x_opt = str.x_opt;

cost = 1;
penalty = 2; % linear penalty
n = 10000; % number of the day simulations

%% Initialize e_l, e_r
pd_e = makedist('Normal',12,sqrt(12));
e = x0 + P*pdf(pd_e,t); % gaussian radiance distribution during the day
x = e;

mu = 0.3; % variance of the gaussian distribution at each moment
e_l = e - mu;
e_u = e + mu;

pd_E = makedist('Normal',0,mu);
E = random(pd_E,n,N); E = E + ones(n,1)*e;

[obj, i_bad, i_good] = calculate_objective;
%% Plot the solutions and data
%% First window
xi = 0:0.1:24;
y_min = 0; y_max = 4;

hMain = figure('SizeChangedFcn', @callback_size_change);
s = hMain.Position;

axMain = axes('Position', [0.05 0.1 0.925 0.75]);
hPanel = uipanel('Title', 'Eingabe', 'Visible', 'off', 'Position', [0 0.9 1 0.1]);           % Eingabepanel
pos = linspace(100,s(3)-50,25);                                             
for k = 1:N;
    hEdit(k) = uicontrol('Style', 'edit', 'Parent', hPanel, 'Visible', 'off', 'String', num2str(x(k)),...
                         'Position', [pos(k) 20 20 20],...
                         'Callback', @update_main_plot);
    hLabel(k) = uicontrol('Style', 'text', 'Parent', hPanel, 'Visible', 'off', 'String', [num2str(k-1) ':00'],...
                          'Position', [pos(k) 40 20 20]);
end
hSlider = uicontrol('Style', 'slider', 'Visible', 'off', 'BackgroundColor', 'white',...
                    'SliderStep', [1/n 100/n], 'Min', 1, 'Max', n, 'Value', 3,...
                    'Position', [0 0 100 20],...
                    'Callback', @update_main_plot);               
hLabel_zeit = uicontrol('Style', 'text', 'Parent', hPanel, 'Visible', 'off', 'String', 'Zeit',...
                        'Position', [0 20 20 20]);
hLabel_wert = uicontrol('Style', 'text', 'Parent', hPanel, 'Visible', 'off', 'String', 'Wert',...
                        'Position', [0 40 20 20]);
hButton = uicontrol('Style', 'pushbutton', 'Parent', hPanel, 'Visible', 'off', 'String', 'START',...
                    'Position', [0 0 0 0],...
                    'Callback', @callback_startbutton);
hButton_xsta = uicontrol('Style', 'pushbutton', 'Visible', 'off', 'String', 'sta',...
                       'Position', [0 0 0 0],...
                       'Callback', @callback_button_xsta);
hButton_xopt = uicontrol('Style', 'pushbutton', 'Visible', 'off', 'String', 'opt',...
                       'Position', [0 0 0 0],...
                       'Callback', @callback_button_xopt);
callback_size_change
%% Callbacks
function callback_startbutton(source, callbackdata)
    x = zeros(1,N);
    for k = 1:N
        object = hEdit(k);
        x(k) = str2double(object.String);
    end
    [obj, i_bad, i_good] = calculate_objective;
    
    update_main_plot
    create_bestworst_plot
    create_general_plot
end
function callback_size_change(source, callbackdata)
    size = hMain.Position(3:4);
    hSlider.Position = [0 0 size(1) 20];
    pos = linspace(50,size(1)-100,25);
    
    hLabel_wert.Position = [0 5 50 20];
    hLabel_zeit.Position = [0 25 50 20];
    hButton.Position = [size(1)-50 5 40 40];
    hButton_xsta.Position = [5 size(2)-85 30 20];
    hButton_xopt.Position = [45 size(2)-85 30 20];
    for k = 1:N;
        object = hEdit(k);
        object.Position = [pos(k) 5 35 20];
        object.Visible = 'on';
        object = hLabel(k);
        object.Position = [pos(k) 25 35 20];
        object.Visible = 'on';
    end
    % Visibility
    hPanel.Visible = 'on';
    hLabel_wert.Visible = 'on';
    hLabel_zeit.Visible = 'on';
    hSlider.Visible = 'on';
    hButton.Visible = 'on';
    hButton_xsta.Visible = 'on';
    hButton_xopt.Visible = 'on';
end
function update_main_plot(source, callbackdata)
    x = zeros(1,N);
    for k = 1:N
        object = hEdit(k);
        x(k) = str2double(object.String);
    end
    [obj, ~, ~] = calculate_objective;
    
    i = round(hSlider.Value);
    yi = pchip(t, E(i,:), xi);      % wether realization (interpolation)
    
    plot(xi, yi, 'k'), hold on      % i-th simulation
    plot(t, x, 'r*',...             % schedule
         [t; t], [e_l; e_u], 'k+-.')% uncertainty intervals
    xlim([0 24])
    ylim([y_min y_max])
    v = axis;                       % for the positioning of the text boxes
    text(0.05*v(2),0.98*v(4), ['Gewinn = ', num2str(obj(i))])
    legend('Sonnenstrahlung waehrend des Tages',...
           'deine Vorhersage')
    xlabel('Zeit'), ylabel('Energie, MWh')
    title(['Simulation ' num2str(i) ' von ' num2str(n)])
    axMain.XLim = [-0.5 24.5];
    axMain.XTick = 0:24;
    grid on
    hold off 
end
function callback_button_xsta(source, callbackdata)
    x = e;
    for k = 1:N
        object = hEdit(k);
        object.String = x(k);
    end
    [obj, i_bad, i_good] = calculate_objective;
    
    update_main_plot
end
function callback_button_xopt(source, callbackdata)
    x = x_opt;
    for k = 1:N
        object = hEdit(k);
        object.String = x(k);
    end
    [obj, i_bad, i_good] = calculate_objective;
    
    update_main_plot
end
function create_bestworst_plot
    figure,
    ax1 = axes;
    yi = pchip(t, E(i_good,:), xi); % best wether realization (interpolation)
    plot(xi, yi, 'b')               % best wether realization (plot)
    hold on
    yi = pchip(t, E(i_bad ,:), xi); % worst wether realization (interpolation)
    plot(xi, yi, 'r')               % worst wether realization (plot)
    plot(t, x, 'k*',...             % schedule
         [t; t], [e_l; e_u], 'k+-.')% uncertainty intervals

    legend('Tag mit dem groessten Gewinn',...
           'Tag mit dem kleinsten Gewinn',...
           'deine Vorhersage')
    xlabel('Zeit'), ylabel('Energie, MWh')
    title('Tage mit dem kleinsten und dem groessten Gewinn')
    % size
    ylim([y_min y_max])
    v = axis; % for the positioning of the text boxes
    % info-box at the top left corner
    text(0.05*v(2),0.98*v(4), ['max Gewinn = ', num2str(obj(i_good))], 'Color','blue')
    text(0.05*v(2),0.95*v(4), ['min Gewinn = ', num2str(obj(i_bad))], 'Color','red')
    text(0.05*v(2),0.90*v(4), ['Preis = ', num2str(cost)])
    text(0.05*v(2),0.87*v(4), ['Strafe = ', num2str(penalty)])
    ax1.XLim = [-0.5 24.5];
    ax1.XTick = 0:24;
    grid on
    hold off
end
function create_general_plot
    figure,
    ax2 = axes;
    % plot
    plot(t, x, 'r*')                               % schedule
    hold on
    plot(t, E(1:40:end,:), 'ok', 'MarkerSize', 1)  % wether realizations
    plot([t; t], [e_l; e_u], 'k+-.')               % uncertainty intervals
    % size
    ylim([y_min y_max])
    v = axis; % for the positioning of the text boxes
    % info-box at the top left corner
    legend('deine Vorhersage',...
           'simulierte Werte')
    xlabel('Zeit'), ylabel('Energie, MWh')
    title('Simulierte Werte und die Vorhersage')
    text(0.05*v(2),0.98*v(4), ['durchschnittlicher Gewinn = ', num2str(mean(obj))])
    ax2.XLim = [-0.5 24.5];
    ax2.XTick = 0:24;
    grid on
    hold off
end
function [obj, i_bad, i_good] = calculate_objective
    X = ones(n,1)*x;
    obj = sum(cost*E - max(zeros(n,N), X - E)*penalty - max(zeros(n,N), E - X)*cost, 2);
    i_bad = find(obj == min(obj), 1, 'first');
    i_good = find(obj == max(obj), 1, 'first');
end

end