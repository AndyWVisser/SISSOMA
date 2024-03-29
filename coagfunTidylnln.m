function sim = coagfunTidylnln(a,alpha,epsilon,Ptotal,Tmax,seasonal)

arguments
    a double = 2.0;         % selfsimilarity parameter; typically between 1.8 and 2.0
    alpha double = 0.3      % stickiness; range 0 to 1
    epsilon double = 1E-6;  % turbulent dissipation rate [m3/s2]
    Ptotal double = 1E6;    % total productivity; typically 1E6 [µg m-2 day-1] (1 gC m-2 day-1)  
    Tmax double = 5*365;     % period of simulation [days]
    seasonal logical = 1;   % seasonal (1) or constant (0) production
end

disp(['a = ',num2str(a,2), ' alpha = ',num2str(alpha,2), ' epsilon = ',num2str(epsilon,2),...
    ' Ptotal = ',num2str(Ptotal,2), ' Tmax = ', num2str(Tmax,5)])

p.a = a; % self-similarity parameter
p.alpha = alpha; % stickiness
p.epsilon = epsilon; % [m^2 s^-3] % energy dissipation rate 
epsilon_crit = 1E-6; % [m^2 s^-3] % energy dissipation rate
p.Ptotal = Ptotal; % [µg C m^2 day^-1] % 
p.seasonal = seasonal;

Frate = 0.5; p.Frate = Frate; % maximum fragmentation rate [day-1] for aggregates > r_crit
H = 20;     p.H = H;        %[m] depth of mixed layer

rMax = 1E6; p.rMax = rMax;  % [µm] maximum radius
ro = 1 ;  p.ro = ro;  %[µm] min radius
r_crit = 1E4; %[µm]

rho_sw = 1027;  p.rho_sw = rho_sw; % density of seawater [kg m^-3]
rhoo = 1;   % minimum aggregate density at r = ro [kg m^-3]
rhom = 200; % maximum aggregate density at r = ro [kg m^-3]

nu = 1E-6;      % [m^2 s^-1] kinematic viscosity of seawater
mu = nu*rho_sw; % [kg m^-1 s^-1] dynamic viscosity of seawater
kb = 1.38065E-23; %Boltzmann constant [m^2 kg s^-2 K^-1]

rfactor = 0.02; % 0.1; p.rfactor = rfactor; % [day^-1] remineralization rate 
MtoC = 2.5; % dry mass to dry mass carbon ratio [µg µgC-1]


%% grid and combination variables
Nr = 30; p.Nr = Nr; %number of size bins
Nd = 10; p.Nd = Nd; %number of density bins

delta = exp(log(rMax/ro)/(Nr-1));  p.delta = delta;   % logrithmic size interval
lambda = exp(log(rhom/rhoo)/(Nd-1)); p.lambda = lambda;   % logrithmic density interval
% drho = 0.1*rho_sw/(Nd-1);          p.drho = drho;     % density interval
q = delta^(p.a-3);
% zw = rho_sw/drho;

L = Nr*Nd;      % number of bins: b index [0, 1, ..., L-1]
K = (L+1)*L/2;  % number of combos: k index [0, 1, ..., K-1]
k = (0:K-1)';
z = (2*L + 1 - sqrt((2*L + 1).^2 - 8*k))/2;%(2*L + 1 - sqrt((2*L + 1).^2 - 8*k))/2;
bi = floor(z);
bj =   k - bi*L + bi.*(bi - 1)/2 + bi;%k - bi*L + bi.*(bi-1)/2 + L (4.13)
xi = floor(bi/Nd); zi = bi - xi*Nd; % Does this make sense? M is nD and N is nR
xj = floor(bj/Nd); zj = bj - xj*Nd;
x = 0:Nr-1; z = 0:Nd-1;

[x_mesh,z_mesh] = meshgrid(x,z);

% functions relating ordinates to values
pip = 4*pi/3;
r_fun = @(x) ro*delta.^x;      % [µm] aggregate radius
v_fun = @(x) pip*r_fun(x).^3;  % [µm^3] aggregate volume
d_fun = @(x,z) rhoo*lambda.^z.*q.^x;    % [kg/m^3] aggregate excess density
p_fun = @(x) delta.^((a-3)*x); % dry mass fraction (phi)
m_fun = @(x,z) v_fun(x).*(d_fun(x,z) + rho_sw)*1E-18*1E9; % [µg] aggregate total mass
md_fun = @(x,z) v_fun(x).*(d_fun(x,z) + p_fun(x)*rho_sw)*1E-18*1E9; % [µg] aggregate dry mass
e_fun = @(x,z) v_fun(x).*d_fun(x,z)*1E-18*1E9; % [µg] aggregate excess mass
w_fun = @(x,z) 9.8*d_fun(x,z).*((r_fun(x)*1E-6).^2)*2./(9*nu*rho_sw)*24*3600 ; % [m/s]*24*3600 [m d^-1] sinking velocity 

r_mean = zeros(size(x));
v_mean = zeros(size(x));
p_mean = zeros(size(x));
d_mean = zeros(size(x_mesh));
m_mean = zeros(size(x_mesh));
md_mean = zeros(size(x_mesh));
e_mean = zeros(size(x_mesh));
w_mean = zeros(size(x_mesh));
for i = 1:Nr
    r_mean(i) = integral(r_fun,i-1,i); % [µm]
    v_mean(i) = integral(v_fun,i-1,i); % [µm]
    p_mean(i) = integral(p_fun,i-1,i); % [µm]
    for j = 1:Nd
        d_mean(j,i) = integral2(d_fun,i-1,i,j-1,j); % [kg/m^3]
        m_mean(j,i) = integral2(m_fun,i-1,i,j-1,j); % [µg]
        md_mean(j,i) = integral2(md_fun,i-1,i,j-1,j); % [µg]
        e_mean(j,i) = integral2(e_fun,i-1,i,j-1,j); % [µg]
        w_mean(j,i) = integral2(w_fun,i-1,i,j-1,j); % [m/day]
    end
end


%% transformations
% xzb = @(b) [floor(b/(Nd)), b - floor(b/(Nd))*(Nd)]; % bin number into x-z coordinates
bxz = @(x,z) x*(Nd) + z; %gives bin number in a vector

% coordinates of the daughter particles for every combination
xioj = log(delta.^(a*xi) + delta.^(a*xj))./log(delta^a);
zioj = log(lambda.^zi.*delta.^(a*(xi-xioj)) + lambda.^zj.*delta.^(a*(xj-xioj)))/log(lambda);

% Target bins
% Dividing mass between four target bins
x300 = floor(xioj); 
z300 = floor(zioj); 
b300 = bxz(x300,z300);    %defining all four target bins (bin number)
b300(b300<0) = 0;
b310 = bxz(x300+1,z300);
b301 = bxz(x300, z300+1);
b311 = bxz(x300+1,z300+1);
dx1 = xioj - x300;    %dividing mass between x and z ordinates
dx0 = 1 - dx1;
dz1 = zioj - z300;
dz0 = 1 - dz1;
f00 = dx0.*dz0; %determining the fraction going into each bin
f10 = dx1.*dz0;
f01 = dx0.*dz1;
f11 = dx1.*dz1;
koutx = x300+1>Nr-1;    % boundary condition: material exiting domain
koutz = z300+1>Nd-1;
f10(koutx) = 0; f11(koutx) = 0; b310(koutx) = b300(koutx); b311(koutx) = b300(koutx);
f01(koutz) = 0; f11(koutz) = 0; b301(koutz) = b300(koutz); b311(koutz) = b300(koutz);


%% environmental variables
T = 281; %temperature


%% sinking speed
w_tmp = w_mean /(24*3600);        % [m/s]
r_m = ones(size(z'))*r_mean*1E-6; %[m]
Re = 2.*r_m.*w_tmp./nu;
Cd = 24./Re + 6./(1+Re.^0.5) + 0.4;
w_it = sqrt((8*r_m.*9.81.*d_mean)./(3*rho_sw.*Cd)); % [m/s]
tick = 0;
while max(w_tmp./w_it,[],'all')>1 %iterate until converged (within 0.5-2 times the previous)
    w_tmp =w_it;
    Re = 2.*r_m.*w_tmp./nu;
    Cd = 24./Re + 6./(1+Re.^0.5) + 0.4;
    w_it = sqrt((8*r_m.*9.81.*d_mean)./(3*rho_sw.*Cd)); % [m/s]
    tick=tick+1;
end
w = w_it*24*3600; % [m / d]


%% coagulation kernels
beta_b = zeros(size(xi));
beta_s = zeros(size(xi));
beta_d = zeros(size(xi));
beta   = zeros(size(xi));
for i = 1:K
    ri = r_mean(xi(i)+1); wi = w_it(xi(i)+1); %[m/s]
    rj = r_mean(xj(i)+1); wj = w_it(xj(i)+1); %[m/s]
    beta_b(i) = (2*kb*T)./(3*mu)*(ri+rj)^2./(ri*rj);  % [m^3/s], Brownian motion  % 
    pp = ri/rj;
    beta_s(i) = 9.8*(pp.^2./(1+2*pp.^2)).*(epsilon/nu)^0.5.*(1E-6*(ri+rj))^3;  %[m^3/s], Turbulent shear
    beta_d(i) = 0.5*pi*abs(wi - wj)*(1E-6*ri).^2; % [m^3/s], Differental settling
    beta(i) = (beta_b(i) + beta_s(i) + beta_d(i))*3600*24; %[m^3 / day ] Total encounter rate
end


%% Porosity estimate
% phi = delta.^((a-3)*x_mesh); % phi == 1 - porosity (i.e. dry mass volume fraction)


%% Fragmentation
% f_fun = @(x) (ro*delta.^x/1E4).^3; % propostional to volume
% pfrag = f_fun(x_mesh).*(1-phi)*(epsilon*1E6);

phi = p_mean.*ones(10,30);

f_fun = @(x) Frate*(r_mean/r_crit).^3; % propostional to volume
pfrag = f_fun(x_mesh).*(1-phi)*(epsilon/epsilon_crit);


%% Remineralization
remin = rfactor/log(lambda)*z_mesh;

% maxd_mean = max(max(d_mean));

% Test
% remin1 = rfactor*((1-q)/(3-a))*(q.^x_mesh).*(1+z_mesh);
% remin2 = rfactor/log(lambda) * (d_mean/maxd_mean);
% remin3 = rfactor*(zw+z_mesh);

%% Interactions
M = zeros(L,1); % [µgC/ m^3]


%% Productivity
prod = zeros(Nd,Nr); 
prod(1:end,1) = MtoC*Ptotal/Nd/H;% production of POM [µg / m^3 / day] (estimate from pimary prodiction)

% mdry = v_mean .* (d_mean * 1E-9 + phi * rho_sw * 1E-9); % [µg] of a singe aggregate
mdry = md_mean;


%% Time dependent solution 
tic
options = odeset('NonNegative',1:length(M(:)));
if seasonal
    %disp('seasonal')
    [t,Mo] = ode15s(@interaxseason, [0:Tmax], [M(:) ],options,mdry,bi,bj,Nr,Nd,b300,b301,b310,b311,f00,f01,f10,f11,alpha,beta,w,H,prod(:),remin(:),rfactor,pfrag(:));
else
    %disp('non-seasonal')
    [t,Mo] = ode15s(@interax, [0 Tmax], [M(:) ],options,mdry,bi,bj,Nr,Nd,b300,b301,b310,b311,f00,f01,f10,f11,alpha,beta,w,H,prod(:),remin(:),rfactor,pfrag(:));
end

runtime = toc  

RMSE = 0;
for i = 1:length(t)-1
    RMSE(i) = rms((Mo(i+1,:)-Mo(i,:)));
end
%% 

% dry mass per volume in last time step in each size-density bin [µgC/m3]
Mdry = reshape(Mo(end,:),Nd,Nr);

% flux out of the H=50m mixed layer [µgC/m2/day]
Flux = Mdry.*w/MtoC;  % flux in µg C m^2 / day

% Size integrated flux [µgC/m2/day]
BFlux = sum(Flux,1);

% number of particles per volume [#/m3]
N = Mdry./mdry;

figure(1); clf;
subplot(2,2,1); imagesc(x,z,Mdry); colorbar;  title('dry mass'); axis xy
subplot(2,2,2); imagesc(x,z,Flux); colorbar; title('flux'); axis xy
subplot(2,2,4); semilogy(t(2:end),RMSE,'.'); title('root mean error');
subplot(2,2,3); plot(x+.5,BFlux,'o'); title('size integrated flux');

sim.prod = prod;
sim.Mdry = Mdry;
sim.Flux = Flux;
sim.RMSE = rms(Mo(end,:)-Mo(end-1,:));
sim.w = w;
sim.N = N;
sim.mTot = m_mean;
sim.phi = phi;
sim.frag = pfrag;
sim.remin = remin;
sim.d = d_mean;
sim.r = r_mean;
sim.md = md_mean;
sim.v = v_mean;
sim.pm = p_mean;
sim.t = t;
sim.Mo = Mo;
sim.p = p;
sim.x = x;
sim.z = z;


[dMdto,dMsink,dMremin,dMfrag] = interax(t,Mdry(:),mdry,bi,bj,Nr,Nd,b300,b301,b310,b311,f00,f01,f10,f11,alpha,beta,w,H,prod(:),remin(:),rfactor,pfrag(:));

sim.dMdto = reshape(dMdto,Nd,Nr);
sim.dMsink = reshape(dMsink,Nd,Nr);
sim.dMremin = reshape(dMremin,Nd,Nr);
sim.dMfrag = reshape(dMfrag,Nd,Nr);

figure(2); clf;
subplot(3,2,1); imagesc(x,z,reshape(dMdto,Nd,Nr)); colorbar; title('dMdt'); axis xy;
subplot(3,2,2); imagesc(x,z,reshape(dMsink,Nd,Nr)); colorbar; title('sinking'); axis xy;
subplot(3,2,4); imagesc(x,z,reshape(dMremin,Nd,Nr)); colorbar; title('remineralization'); axis xy;
subplot(3,2,3); imagesc(x,z,reshape(dMfrag,Nd,Nr)); colorbar; title('fragmentation'); axis xy;
subplot(3,2,5); imagesc(x,z,prod); colorbar; title('production'); axis xy;
end





