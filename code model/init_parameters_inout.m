function [T, P, cost, penalty, penalty_grad, epsilon, C, SOC_0, t, mu, sigma, lambda, penalty_hess] = init_parameters_inout

% T = 13;
%T = 1 + 24*1;       % one hour schedule (works for RO and start_SO_closed, for the latter only with long running time)
% T = 3;              % for start_SO_closed with an acceptable runtime
% T = 1440;           % every-minute schedule (required for start_SO_discretization and start_SO_discretization_battery when using PVdata2)
 T = 96;              % 15min-schedule (required for start_SO_discretization and start_SO_discretization_battery when using sample_normal_independent)
         
P = 3.8;             % nominal power of the PV element
epsilon = 0.05;
cost = ones(1,T)*5;   % cost(j) is a price for an energy unit during j's hour
load('realistic pricing')
% cost(18*4:22*4) = 2*cost(18*4:22*4);

penalty = @(x) x.^2;  % quadratic penalty
penalty_grad = @(x) 2*x; % derivative of the penalty function
penalty_hess = @(x) 2; % Hessian of the penalty function
% penalty = @(x) 10*x; % linear penalty
% penalty_grad = @(x) 10*ones(size(x)); % derivative of the penalty function
% penalty_hess = @(x) zeros(size(x)); % Hessian of the penalty function
C = 2.6;              % battery capacity
SOC_0 = 0.25;         % state of charge (SOC) at day break
% SOC_0 = 0.7;
t = linspace(0,24,T);

% pd = makedist('Normal');
% mu = P*(pdf(pd,(t-12)/sqrt(12)));
% sigma = 0.2*ones(1,T);

% Estimated parameters for normal distribution
mu = 1e-3*[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0716994390681004,1.75638989964158,30.1538241236559,87.5189526397849,157.329069449462,240.279084712545,397.984206810036,589.491841487455,774.452829928315,935.517520071685,1078.95311308244,1222.16588512545,1348.3307483871,1465.00860107527,1608.59112419355,1768.06379193548,1851.23593530466,1936.28961845878,2017.40100448029,2082.39888781362,2202.70638691756,2227.48831612903,2262.70253888889,2226.58467275986,2260.25099623656,2294.42242759857,2297.45506379928,2296.49896827957,2274.33150896057,2273.50836648746,2239.56278853047,2167.31725071685,2070.91292258065,1992.72602168459,1940.11473924731,1891.58183297491,1809.80839767025,1686.85554265233,1540.47001379928,1445.89327329749,1329.36601057348,1217.15813691756,1033.41883225806,883.583432616487,732.163161917563,566.027376702509,383.974869964158,260.075257849462,170.406198324373,86.9768640053763,14.3391170716846,3.26036262724014,0.262570605734767,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
mu = mu(floor(linspace(1,length(mu),T)));
sigma = 1e-3*[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.301841439722851,4.77171908278249,41.1954667230687,101.586231411123,159.1036038286,221.183932203628,268.1978031226,314.903965665389,365.575742259369,435.485644076735,494.351486385903,545.151966253382,593.445223194853,639.352201127092,675.459578140996,706.089060838083,739.898572488374,794.186238691072,857.826313801539,858.97252755491,840.726054531225,835.936265554767,840.430597638838,881.245821698738,888.322252299938,893.738109936579,904.977402018073,890.213724088489,875.153514178022,884.602131925548,869.041859719286,828.547638499578,846.101782522385,849.427279623307,808.944890221886,761.884565024479,730.233824262124,713.236190715957,673.246742480306,617.764147572443,588.335070889555,529.118623534293,479.446690702823,411.395784140893,377.620853338791,342.859111663233,301.012211249166,222.879167925889,170.639874604181,102.508542513154,29.9221279065132,10.0532672876217,1.15524268748841,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
sigma = sigma(floor(linspace(1,length(sigma),T)));

% Estimated parameters for weibull distribution
% shape = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.8443078447933,0.708217014772389,0.833049673062412,0.737391486508102,0.526191784713442,0.87131752783714,1.45875443169578,1.94061886516121,2.33925262181793,2.40404840539718,2.48101246371973,2.5304906845987,2.55077543950141,2.6354461043423,2.85898792140101,3.07615658673553,3.03231367535781,2.93233273413678,2.80844251022367,2.91775474924801,3.28237274022606,3.39515374629108,3.45343790356611,3.12896584816938,3.14246447741732,3.18784687285897,3.11778904849928,3.20782848494289,3.25660414629743,3.23474361751589,3.26421684507622,3.29474231874623,2.95426085204002,2.78059325644164,2.88768928015204,3.04278511931806,2.99296828175192,2.77821548750608,2.66416226502895,2.70620087633985,2.52037493272543,2.62418606961156,2.42607492108726,2.34825881227957,2.03177874297603,1.63907857392463,1.05374700068333,0.767424726790466,0.689458153055244,0.845287652840669,0.692747682405154,0.622622954776224,0.718635407305568,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
% scale = [1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,0.429027703327558,3.79792441171697,44.8441135382855,101.121517015507,111.058794318293,237.016079945169,445.277323024775,671.107572916055,884.892295436995,1070.18446214169,1233.70283114668,1392.07734072158,1533.63307968473,1669.47782087834,1836.17324732368,2010.36812464203,2103.33069196264,2203.82566554949,2300.37139881932,2372.73480007181,2497.16751364429,2525.01078585557,2562.03762468268,2533.4842084137,2570.12355377665,2604.94769216863,2606.69907803909,2608.30439956599,2582.16795136131,2581.54055801467,2545.13986906805,2460.39998060122,2358.22197607585,2273.93277695913,2212.21266691118,2153.78075889733,2056.58765283767,1919.05431043029,1758.17695695131,1645.18543927697,1510.12106276439,1384.12416197627,1181.5087238052,1006.45469812403,833.863257175453,639.249915039398,398.356958058637,251.890024227151,163.917292866299,112.473858091869,18.4448104302872,5.63494264557899,1.17430911499071,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06];
% mu = 1e-3*gamma(1+1./shape) .* scale;
% sigma = 1e-6*(gamma(1+2./shape) - gamma(1+1./shape).^2) .* (scale.^2);

lambda = 1.645; % lambda*sigma is the width of the 90%-confidence intervals in RO (for normal distribution)
% lambda = 0.675; % lambda*sigma is the width of the 50%-confidence intervals in RO (for normal distribution)