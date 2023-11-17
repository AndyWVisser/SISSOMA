rmins = [1 1];
rmaxs = [1E5 1E6];

for i = 1:length(rmins)

    a = 2;           % selfsimilarity parameter; typically between 1.8 and 2.0
    alpha = 0.4;     % stickiness; range 0 to 1
    epsilon = 1E-6;  % turbulent dissipation rate [m3/s2]
    Ptotal = 1E6;    % total productivity; typically 1E6 [µg m-2 day-1] (1 gC m-2 day-1)
    Rrate = 0.1;     % remineralization rate [day-1]
    Frate = 500;     % maximum fragmentation rate [day-1] for aggregates > 1 m
    Tmax = 5*365;    % period of simulation [days]
    seasonal = 1;   % seasonal (true) or constant (false) production

    lowlim = 15;
    uplim = 20;

    rMax = rmaxs(i);
    ro = rmins(i);

    T_input = 10;

    [r_mean, d_mean, v_mean, m_mean, w, pfrag, remin, phi] = compareFUN(a,alpha,epsilon,Ptotal,Rrate,Frate,Tmax,seasonal,T_input,ro,rMax);

    r_mean2(i,:) = r_mean;
    d_mean2(i,:,:) = d_mean;
    m_mean2(i,:,:) = m_mean;
    w2(i,:,:) = w;
    pfrag2(i,:,:) = pfrag;
    remin2(i,:,:) = remin;
    phi2(i,:,:) = phi;


end

%%
col = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.3010 0.7450 0.1330], [0.1010 0.2450 0.9330], [0.7010 0.1450 0.9330], [0.3010 0.7450 0.1330] };

param = {d_mean2, w2, pfrag2, remin2, phi2};
tit = {'excess density', 'sinking speed', 'fragmentation', 'remineralization', 'phi'};

% 
for j = 1:length(param)
    
    figure
    for i = 1:10
        loglog(squeeze(r_mean2(1,:,:)), squeeze(param{j}(1,i,:)), '.-', 'Color', col{i}, 'Linewidth', 0.5)
        hold on
        loglog(squeeze(r_mean2(2,:,:)), squeeze(param{j}(2,i,:)), '--', 'Color', col{i}, 'Linewidth', 1)

    end

    title(tit{j})
    hold off
    xlabel('radius [\mum]')
    ylabel('rate [1/day]')
end


%%

interp1(squeeze(r_mean2(2,:,:)), squeeze(param{4}(2,i,:)), 10)
interp1(squeeze(r_mean2(2,:,:)), squeeze(param{4}(1,i,:)), 10)

%%
Nr = 1:30;
Nd = 1:10;

figure
imagesc(remin2, Nd, squeeze(d_mean2(2,:,:))./squeeze(d_mean2(1,:,:))); colorbar; title('excess density'); axis xy

figure
imagesc(Nr, Nd, squeeze(m_mean2(2,:,:))./squeeze(m_mean2(1,:,:))); colorbar; title('total mass'); axis xy

figure
imagesc(Nr, Nd, squeeze(w2(2,:,:))./squeeze(w2(1,:,:))); colorbar; title('sinking speed'); axis xy

figure
imagesc(Nr, Nd, squeeze(pfrag2(2,:,:))./squeeze(pfrag2(1,:,:))); colorbar; title('fragmentation rate'); axis xy

figure
imagesc(Nr, Nd, squeeze(remin2(2,:,:))./squeeze(remin2(1,:,:))); colorbar; title('remineralization rate'); axis xy

%%
figure
semilogx(squeeze(r_mean2(2,:,:)), ones(1,30), 'o', 'Linewidth', 1.2); title('size grid');
hold on
semilogx(squeeze(r_mean2(1,:,:)), ones(1,30), 'x', 'Linewidth', 1.2); title('size grid');

%%
figure
imagesc(Nr, Nd, squeeze(phi2(2,:,:))./squeeze(phi2(1,:,:))); colorbar; title('phi'); axis xy



%%


function [r_mean, d_mean, v_mean, m_mean, w, pfrag, remin, phi] = compareFUN(a,alpha,epsilon,Ptotal,Rrate,Frate,Tmax,seasonal,T_input, r0, r_max)

arguments
    a double = 2;           % selfsimilarity parameter; typically between 1.8 and 2.0
    alpha double = 0.4      % stickiness; range 0 to 1
    epsilon double = 1E-6;  % turbulent dissipation rate [m3/s2]
    Ptotal double = 1E6;    % total productivity; typically 1E6 [µg m-2 day-1] (1 gC m-2 day-1)
    % Ptotal double = 0.03*1E3*50;    % total productivity; typically 1E6 [µg m-2 day-1] (1 gC m-2 day-1)
    Rrate double = 0.1;     % remineralization rate [day-1]
    Frate double = 500;     % maximum fragmentation rate [day-1] for aggregates > 1 m
    Tmax double = 5*365;    % period of simulation [days]
    seasonal logical = 0;   % seasonal (true) or constant (false) production
    % temp_depend_remin logical = 1; % constant remineralization rate (false) or temperature dependent (true)
    T_input double = 10;
    r0 double =[];
    r_max double =[];
end

p.a = a; % self-similarity parameter
p.alpha = alpha; % stickiness
p.epsilon = epsilon; % [m^2 s^-3] % energy dissipation rate
p.Ptotal = Ptotal; % [µg C m^2 day^-1] %
p.seasonal = seasonal;
H = 50;     p.H = H;        % [m] depth of mixed layer
strata = 1; p.strata = strata; % not used

% rMax = 1E6; p.rMax = rMax;  % [µm] maximum radius
% ro = 1 ;  p.ro = ro;  % [µm] min radius
rMax = r_max; p.rMax = rMax;  % [µm] maximum radius
ro = r0;  p.ro = ro;  % [µm] min radius


% rho_sw = 1027;  p.rho_sw = rho_sw; % density of seawater [kg m^-3]
% nu = 1E-6;      % [m^2 s^-1] kinematic viscosity of seawater
% mu = nu*rho_sw; % [kg m^-1 s^-1] dynamic viscosity of seawater

if seasonal == 0 && length(T_input) > 1
    T_input = mean(T_input);
end

rho_sw = zeros(length(T_input),1);
nu = zeros(length(T_input),1);
mu = zeros(length(T_input),1);

for j = 1:length(T_input)
    rho_sw(j) =  SW_Density(T_input(j), 'C', 35, 'ppt', 3, 'bar');
    nu(j) = SW_Kviscosity(T_input(j), 'C', 35, 'ppt');
    mu(j) = nu(j)*rho_sw(j);
end

rho_sw = mean(rho_sw);
nu = mean(nu);
mu = mean(mu);

kb = 1.38065E-23; % Boltzmann constant [m^2 kg s^-2 K^-1]
p.Rrate = Rrate; % [day^-1] remineralization rate
p.Frate = Frate; % [day^-1] maximum fragmentation rate
epsilon_ref = 1E6; p.epsilon_ref = epsilon_ref; % ?????? [s2/m3]

%% grid and combination variables
Nr = 30; p.Nr = Nr; %number of size bins
Nd = 10; p.Nd = Nd; %number of density bins

delta = exp(log(rMax/ro)/(Nr-1));  p.delta = delta;   % logrithmic size interval
drho = 0.2*rho_sw/(Nd-1);          p.drho = drho;     % density interval
q = delta^(p.a-3);

L = Nr*Nd;      % number of bins: b index [0, 1, ..., L-1]
K = (L+1)*L/2;  % number of combos: k index [0, 1, ..., K-1]
k = [0:K-1]';
b = [0:L-1]';
z = (2*L + 1 - sqrt((2*L + 1).^2 - 8*k))/2;
bi = floor(z);
bj =   k - bi*L + bi.*(bi - 1)/2 + bi; %k - bi*L + bi.*(bi-1)/2 + L (4.13)
xi = floor(bi/Nd); zi = bi - xi*Nd; % Does this make sense? M is nD and N is nR
xj = floor(bj/Nd); zj = bj - xj*Nd;
x = [0:Nr-1]; z = [0:Nd-1];

% functions relating ordinates to values
pip = 4*pi/3;
r_fun = @(x) ro*delta.^x;      % [µm] aggregate radius
d_fun = @(x,z) drho*z.*q.^x;   % [kg/m^3] aggregate excess density
v_fun = @(x,z) pip*r_fun(x).^3;% [µm^3] aggregate volume
m_fun = @(x,z) v_fun(x).*(d_fun(x,z) + rho_sw)*1E-18*1E9; % [µg] aggregate total mass
e_fun = @(x,z) v_fun(x).*d_fun(x,z)*1E-18*1E9; % [µg] aggregate excess mass
w_fun = @(x,z) 9.8*d_fun(x,z).*((r_fun(x)*1E-6).^2)*2./(9*nu*rho_sw)*24*3600 ; % [m/s]*24*3600 [m d^-1] sinking velocity

[x_mesh,z_mesh] = meshgrid(x,z);

r_mean = zeros(size(x));
d_mean = zeros(size(x_mesh));
v_mean = zeros(size(x_mesh));
m_mean = zeros(size(x_mesh));
e_mean = zeros(size(x_mesh));
w_mean = zeros(size(x_mesh));
for i = 1:Nr
    r_mean(i) = integral(r_fun,i-2,i-1); % [µm]
    for j = 1:Nd
        d_mean(j,i) = integral2(d_fun,i-1,i,j-1,j); % [kg/m^3]
        v_mean(j,i) = integral2(v_fun,i-1,i,j-1,j); % [kg/m^3]
        m_mean(j,i) = integral2(m_fun,i-1,i,j-1,j); % [µg]
        e_mean(j,i) = integral2(e_fun,i-1,i,j-1,j); % [µg]
        w_mean(j,i) = integral2(w_fun,i-1,i,j-1,j); % [m/day]
    end
end


%% transformations
xzb = @(b) [floor(b/(Nd)), b - floor(b/(Nd))*(Nd)]; % bin number into x-z coordinates
logd = @(x) log(x)./log(delta); %used for finding daughter particles
zeta = @(xi,xj) (1 + delta.^(a*(xi - xj)));
bxz = @(x,z) x*(Nd) + z; %gives bin number in a vector

%coordinates of the daughter particles for every combination
xioj = xi + logd(zeta(xj,xi))/a;
zioj = zi./zeta(xj,xi) + zj./zeta(xi,xj);

% Target bins
% Dividing mass between four target bins
x300 = floor(xioj);
z300 = floor(zioj);
b300 = bxz(x300,z300);    %defining all four target bins (bin number)
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
while max(w_tmp./w_it,[],'all')>1 % iterate until converged (within 0.5-2 times the previous)
    w_tmp = w_it;
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
    beta_b(i) = (2*kb*T)./(3*mu)*(ri+rj)^2./(ri*rj);  % [m^3/s], Brownian motion
    pp = ri/rj;
    beta_s(i) = 9.8*(pp.^2./(1+2*pp.^2)).*(epsilon/nu)^0.5.*(1E-6*(ri+rj))^3;  %[m^3/s], Turbulent shear
    beta_d(i) = 0.5*pi*abs(wi - wj)*(1E-6*ri).^2; % [m^3/s], Differental settling
    beta(i) = (beta_b(i) + beta_s(i) + beta_d(i))*3600*24; %[m^3 / day ] Total encounter rate
end


%% Porosity estimate
phi = delta.^((a-3)*x_mesh); % phi == 1 - porosity (i.e. dry mass volume fraction)


%% Fragmentation rate [1/day]

% Old 
% f_fun = @(x,z) (delta.^x/(delta^(Nr-1))); % proportional to r
% pfrag = Frate*f_fun(x_mesh,z_mesh).*(1-phi);
% pfrag = pfrag*epsilon*epsilon_ref;




% f_fun = @(x,z) (delta.^x/(delta^(Nr-1))); % proportional to r
% f_fun = @(x,z) delta.^x.*1E-6; % proportional to r
% pfrag = f_fun(x_mesh,z_mesh);



% New
f_fun = @(x) delta.^x.*1E-6; % proportional to r
pfrag = Frate*f_fun(x_mesh).*(1-phi);
pfrag = pfrag*epsilon*epsilon_ref;

%% Remineralization

Rrate_ref10 = 0.07; % remineralisation rate (1/day) (Serra-Pompei (2022)) @10 degrees
Q10 = 2;
Rrate = Rrate_ref10 .* fTemp(Q10, T_input);



% if temp_depend_remin == 1 % Make it temperature-dependent:

% if seasonal == 1
%     % daily tmps from TMs (!!! get them from the Transport matrix)
%     % T_mat = [10 10 linspace(10,17,182) linspace(17,10,182) linspace(10,17,182) linspace(17,10,182) linspace(10,17,182) linspace(17,10,182) linspace(10,17,182) linspace(17,10,182) linspace(10,17,182) linspace(17,10,182) 10 10 10 ];
%     % T_mat = [10 10 linspace(10,20,182) linspace(20,10,182) linspace(10,20,182) linspace(20,10,182) linspace(10,20,182) linspace(20,10,182) linspace(10,20,182) linspace(20,10,182) linspace(10,20,182) linspace(20,10,182) 10 10 10 ];
%     Q10 = 2;
%     Rrate = Rrate_ref10 .* fTemp(Q10, T_input);
%
% elseif seasonal == 0
%     Q10 = 2;
%     Rrate = Rrate_ref10 .* fTemp(Q10, T_input);
%
% end
%end

remin_grid = ((1-q)/(3-a)).*(q.^x_mesh).*(1+z_mesh);

remin = Rrate.*remin_grid;


%% Interactions
% N = zeros(L,1); % number of particles/m^3
M = zeros(L,1); % [\mug / m^3]

prod = zeros(Nd,Nr);
prod(1:end,1) = Ptotal/Nd/H; % [\mugC / m^3 / day]


mdry = v_mean .* (d_mean * 1E-9 + phi * rho_sw * 1E-9); % [µm] of a singe aggregate




end






