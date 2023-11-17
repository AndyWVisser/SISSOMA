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
v_fun = @(x,z) pip*r_fun(x).^3;% [µm^3] aggregate volume

%%
iT = 1:length(T_input);
% for i = 1:length(T_input)
    d_fun = @(x,z) drho.*z.*q.^x;   % [kg/m^3] aggregate excess density
%     m_fun = @(iT,x,z) v_fun(x).*(d_fun(iT,x,z) + rho_sw(iT))*1E-18*1E9; % [µg] aggregate total mass
%     e_fun = @(iT,x,z) v_fun(x).*d_fun(iT,x,z)*1E-18*1E9; % [µg] aggregate excess mass
%     w_fun = @(iT,x,z) 9.8*d_fun(iT,x,y).*((r_fun(x)*1E-6).^2)*2./(9*nu(iT)*rho_sw(iT))*24*3600 ; % [m/s]*24*3600 [m d^-1] sinking velocity
% end

%%

[x_mesh,z_mesh] = meshgrid(x,z);

% r_mean = zeros(length(T_input),size(x));
d_mean = zeros(length(T_input),10,30);
% v_mean = zeros(length(T_input),size(x_mesh));
% m_mean = zeros(length(T_input),size(x_mesh));
% e_mean = zeros(length(T_input),size(x_mesh));
% w_mean = zeros(length(T_input),size(x_mesh));
tic
for iT = 1:length(T_input)
    for i = 1:Nr
        r_mean(i) = integral(r_fun,i-2,i-1); % [µm]
        for j = 1:Nd
            d_mean(iT,j,i) = integral2(d_fun,i-1,i,j-1,j); % [kg/m^3]
%             v_mean(j,i) = integral2(v_fun,i-1,i,j-1,j); % [kg/m^3]
%             m_mean(iT,j,i) = integral3(m_fun,iT,iT,i-1,i,j-1,j); % [µg]
% % %             e_mean(iT,j,i) = integral3(e_fun,iT,iT,i-1,i,j-1,j); % [µg]
%             w_mean(iT,j,i) = integral3(w_fun,iT,iT,i-1,i,j-1,j); % [m/day]
        end
    end
end

toc
%%
iT = 1:length(T_input);
iR = 1:Nr;
iD = 1:Nd;

d_mean(iT,j,i) = integral3(d_fun,iT,iT,i-1,i,j-1,j);




%%

    % for iTime = 1:length(T_input)
    %     for i = 1:Nr
    %         r_mean{iTime} = integral(r_fun,i-2,i-1); % [µm]
    %         for j = 1:Nd
    %             d_mean{iTime} = integral2(d_fun{iTime},i-1,i,j-1,j); % [kg/m^3]
    %             v_mean{iTime} = integral2(v_fun{iTime},i-1,i,j-1,j); % [kg/m^3]
    %             m_mean{iTime} = integral2(m_fun{iTime},i-1,i,j-1,j); % [µg]
    %             e_mean{iTime} = integral2(e_fun{iTime},i-1,i,j-1,j); % [µg]
    %             w_mean{iTime} = integral2(w_fun{iTime},i-1,i,j-1,j); % [m/day]
    %         end
    %     end
    % end