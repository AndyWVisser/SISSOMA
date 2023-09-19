
function sim = coagfun(a,alpha,epsilon,Ptotal,simo,Tmax,seasonal)

arguments
    a double            % selfsimilarity parameter; typically between 1.8 and 2.0
    alpha double        % stickiness; range 0 to 1
    epsilon double      % turbulent dissipation [m3/s2]
    Ptotal double       % total productivity; typically 1E5 [µgC m-2 day-1] (0.1 gC m-2 day-1)
    simo struct         % input simulation structure; blank [] for startup
    Tmax double         % period of simulation [days]
    seasonal logical    % seasonal (true) or constant (false) production
end


p.a = a; %fractal dimension
p.alpha = alpha; %stickiness
p.epsilon = epsilon; % [m^2 s^-3] %energy dissipation rate (McCave 1984) (converted from 1E-4 cm^2 s^-3)(1E-8)
p.Ptotal = Ptotal; % [µg C m^2 day^-1] % roughly 10% of net primary production
H = 50;     p.H = H;        %[m] depth of mixed layer
rMax = 1E6; p.rMax = rMax;  %[µm] max radius
rMin = 1 ;  p.rMin = rMin;   %[\µm] min radius
rho_sw = 1.027E-6; % density of seawater [µg µm^-3] (from andy)
nu = 1E-6; % [m^2 s^-1] kinematic viscosity of seawater (from andy)
mu = nu*rho_sw*10^9;% [kg m^-1 s^-1] absolute viscosity (10^9 is a conversion factor for rho to kg/m^3) 
kb = 1.38065E-23; %Boltzmann constant [m^2 kg s^-2 K^-1]
remin = 0.1; p.remin = remin; % [day^-1] remineralization rate 



%% grid and combination variables
nR = 20; %number of size bins
nD = 30; %number of density bins

delta = exp(log(rMax/rMin)/(nR-1));  p.delta = delta;  % size step
deltaRho = 0.8*rho_sw/(nD-1);        p.drho = deltaRho;   % 0.6*rho_sw/nD; %density step
q = delta^(p.a-3); % is this still valid when a is not part of delta?

L = nR*nD; %number of bins: b index [0, 1, ..., L-1]
K = (L+1)*L/2; % number of combos: k index [0, 1, ..., K-1]
k = [0:K-1]';
b = [0:L-1]';
z = (2*L + 1 - sqrt((2*L + 1).^2 - 8*k))/2;%(2*L + 1 - sqrt((2*L + 1).^2 - 8*k))/2; % SHOULD IT BE 2L or L? is this eq 4.12?
bi = floor(z);
bj =   k - bi*L + bi.*(bi - 1)/2 + bi;%k - bi*L + bi.*(bi-1)/2 + L (4.13)
xi = floor(bi/nD); zi = bi - xi*nD; % Does this make sense? M is nD and N is nR
xj = floor(bj/nD); zj = bj - xj*nD;
x = [0:nR-1]; z = [0:nD-1];
[xMesh,zMesh] = meshgrid(x,z);

%% transformations

pip = 4*pi/3;

xz = @(b) [floor(b/(nD)), b - floor(b/(nD))*(nD)]; % bin number into x-z coordinates
logd = @(x) log(x)./log(delta); %used for finding daughter particles
zeta = @(xi,xj) (1 + delta.^(a*(xi - xj)));
pxz = @(x,z) x*(nD) + z; %gives bin number in a vector

r = @(x) rMin*delta.^x; % radius from x ordinate
y = @(x,z) deltaRho*z.*q.^x; % rho-rho_sw from x and z ordinate
mass = @(x,z) pip*(y(x,z) + rho_sw).*(r(x).^3);

w_func = @(x,z) (2.*1E9*y(x,z).*9.81.*(r(x)*1E-6).^2)./(9*mu);% [m/s]   *24*3600; %[m d^-1] sinking velocity NB! from de la rocha and passow 2007, unit double checked

%coordinates of the daughter particles for every combination
xioj = xi + logd(zeta(xj,xi))/a;
zioj = zi./zeta(xj,xi) + zj./zeta(xi,xj);

%% Target bins
% Dividing mass between four target bins
x300 = floor(xioj); x300b = min(x300,nR-1); %lowest x ordinate; 
z300 = floor(zioj); z300b = min(z300,nD-1); %lowest z ordinate
b300 = pxz(x300b,z300b);    %defining all four target bins (bin number)
b310 = pxz(x300b+1,z300b);
b301 = pxz(x300b, z300b+1);
b311 = pxz(x300b+1,z300b+1);
dx1 = xioj - x300b;    %dividing mass between x and z ordinates
dx0 = 1 - dx1;
dz1 = zioj - z300b;
dz0 = 1 - dz1;
f00 = dx0.*dz0; %determining the fraction going into each bin
f10 = dx1.*dz0;
f01 = dx0.*dz1;
f11 = dx1.*dz1;
k01 = f01<0; h10 = b310 > (L-1); h01 = b301 > (L-1); h11 = b311 > (L-1);
k00 = f00<0;
f10(k00) = f10(k00) + f00(k00); f00(k00) = 0; b310(h10) = b300(h10); b301(h01) = L-1;
f11(k01) = f11(k01) + f01(k01); f01(k01) = 0; b311(k01) = b301(k01); b311(h11) = b301(h11);

%% environmental variables
T = 281; %temperature

%% derived properties
%w = @(r,rho) (2.*1E9*(rho-rho_sw).*9.81.*(r*1E-6).^2)./(9*mu);% [m/s]   *24*3600; %[m d^-1] sinking velocity NB! from de la rocha and passow 2007, unit double checked
W = w_func(xMesh,zMesh)*24*3600;%[m d^-1]
    


w_tmp = w_func(xMesh,zMesh);
Re = 2E-6.*r(xMesh).*w_tmp./nu;
fRe = 24./Re + 6./(1+Re.^0.5)+0.4;
w_it = sqrt((8E-6*r(xMesh).*9.81*1E9.*y(xMesh,zMesh))./(3E9*rho_sw.*fRe));
tick = 0;
while max(w_tmp./w_it,[],'all')>1 %iterate until converged (within 0.5-2 times the previous)
    w_tmp =w_it;
    Re = 2E-6.*r(xMesh).*w_it./nu;
    fRe = 24./Re + 6./(1+Re.^0.5)+0.4;
    w_it = sqrt((8E-6*r(xMesh).*9.81*1E9.*y(xMesh,zMesh+.5))./(3E9*rho_sw.*fRe));
    
    
    tick=tick+1;
end
wWhites = w_it*24*3600; %[m/d]



% Vector version for parent particles, used for differential settling
for i = 1:K
    wVeci(i,:) = w_it(zi(i)+1,xi(i)+1);
    wVecj(i,:) = w_it(zj(i)+1,xj(i)+1);
end





%% coagulation kernels
% Brownian motion
beta_b = (2*kb*T)./(3*mu)*(r(xi)+r(xj)).^2./(r(xi).*r(xj));  % [m^3/s], double checked

% Shear
pp = r(xi)./r(xj);
beta_s = 9.8*(pp.^2./(1+2*pp.^2)).*(epsilon/nu)^0.5.*(1E-6*(r(xi)+r(xj))).^3;  %[m^3/s], double checked

% differential settling
beta_d = 0.5*pi*(1E-6*r(xi)).^2.*abs(wVeci-wVecj); % [m^3/s], double checked
 
beta = (beta_b + beta_s + beta_d)*3600*24; %[m^3 d^-1]



%% Fragmentation


pfrag = linspace(0.001,0.5,nR);
pfrag = repmat(pfrag,nD,1);
pfrag = pfrag./(rho_sw+y(xMesh,zMesh))*1E-6;
denfrag = linspace(1,0.5,nD);
for i = 1:nR
    pfrag(:,i) = pfrag(:,i).*denfrag';
end
pfrag = pfrag(:);
frag_div = 0.5;

%% Interactions
N = zeros(nD,nR); %number of particles/m^3 
M = zeros(nD,nR); % [\mug C/ m^3]
m = mass(xMesh,zMesh);
% prod_tot = 1E5; % 0.1 g/m2/d
prod = zeros(size(M));

prod(1:4,1) = Ptotal/9;
prod(5,2:3) = Ptotal/9;
prod(end,3:5) = Ptotal/9;
N = M./m;

%% Time dependent solution 
tic
options = odeset('NonNegative',1:length(M(:)));
if seasonal
    disp('seasonal')
    [t,dM] = ode23(@interactionsDTseason, [0:Tmax], [M(:) ],options,m,xz,bi,bj,nR,nD,q,a,b300,b301,b310,b311,f00,f01,f10,f11,alpha,beta,wWhites,L,H,prod,remin,pfrag,frag_div);
else
    disp('non-seasonal')
    [t,dM] = ode23(@interactionsDT, [0:Tmax], [M(:) ],options,m,xz,bi,bj,nR,nD,q,a,b300,b301,b310,b311,f00,f01,f10,f11,alpha,beta,wWhites,L,H,prod,remin,pfrag,frag_div);
end

runtime = toc  

sim.prod = prod;
sim.M = reshape(dM(end,:),nD,nR);
sim.RMSE = rms(dM(end,:)-dM(end-1,:));
sim.w = wWhites;
sim.N = sim.M./m;
sim.m = m;
sim.d = y(xMesh,zMesh);
sim.r = r(x);
sim.t = t;
sim.dM = dM;
sim.p = p;
sim.x = x;
sim.z = z;
sim.frag = reshape(pfrag,nD,nR);









