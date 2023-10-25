a = 2;           % selfsimilarity parameter; typically between 1.8 and 2.0
alpha = 0.4;     % stickiness; range 0 to 1
epsilon = 1E-6;  % turbulent dissipation rate [m3/s2]
Ptotal = 1E6;    % total productivity; typically 1E6 [µg m-2 day-1] (1 gC m-2 day-1)
Rrate = 0.1;     % remineralization rate [day-1]
Frate = 500;     % maximum fragmentation rate [day-1] for aggregates > 1 m
Tmax = 5*365;    % period of simulation [days]
seasonal = 1;   % seasonal (true) or constant (false) production

lowlim = 10;
uplim = 20;

% if you give a vector of Temperatures --> Temperature dependence (take the mean of the range):
% - rho_sw
% - mu
% - nu
% - coagulation kernels --> check Brownian motion T parameter
% - Rrate

T_input = 10;

% To use the Temperature sesonality --> make sure to follow the productivity pattern!!!
% T0 = 10;
% tt = 1:Tmax;
% T_input = T0*(3-cos(2*pi*tt/365))/3;
% -------------------------------------
% T_input = [lowlim lowlim linspace(lowlim,uplim,182) linspace(uplim,lowlim,182) linspace(lowlim,uplim,182) linspace(uplim,lowlim,182) linspace(lowlim,uplim,182) linspace(uplim,lowlim,182) linspace(lowlim,uplim,182) linspace(uplim,lowlim,182) linspace(lowlim,uplim,182) linspace(uplim,lowlim,182) lowlim lowlim lowlim ];

sim = coagfunDev(a,alpha,epsilon,Ptotal,Rrate,Frate,Tmax,seasonal,T_input);


% !!!!!!!!!!!!!

% - check the interpolation in interaxseason

% - change of the drho for different temperatures

% - Warning: Failure at t=1.824000e+03.  Unable to meet integration tolerances without reducing the step size below the
% smallest value allowed (3.637979e-12) at time t. 
