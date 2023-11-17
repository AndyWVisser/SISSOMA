rMax = 1E5; p.rMax = rMax;  % [µm] maximum radius
ro = 1;  p.ro = ro;  % [µm] min radius

rho_sw = 1027;  p.rho_sw = rho_sw; % density of seawater [kg m^-3]
nu = 1E-6;      % [m^2 s^-1] kinematic viscosity of seawater
mu = nu*rho_sw; % [kg m^-1 s^-1] dynamic viscosity of seawater

a = 2;


Nr = 30; p.Nr = Nr; %number of size bins
Nd = 10; p.Nd = Nd; %number of density bins

delta = exp(log(rMax/ro)/(Nr-1));  p.delta = delta;   % logrithmic size interval
drho = 0.2*rho_sw/(Nd-1);          p.drho = drho;     % density interval
q = delta^(a-3);

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


%%
figure
x = ones(1,30);

plot(r6, x, 'o')
hold on
plot(r5, x,'x')


%%


x6 = w;
x5 = load("x5_w.mat");
x5 = cell2mat(struct2cell(x5));

x = x6./x5;

figure;
imagesc(Nr, Nd, x);
colorbar





%%
% x = 0.0065;
% mugC = x * 12 * 1E6
% model = sum(sim.BFlux)


