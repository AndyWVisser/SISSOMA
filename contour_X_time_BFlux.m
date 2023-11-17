Mdry = zeros(length(sim.t),sim.Nd,sim.Nr);
Flux = zeros(length(sim.t),sim.Nd,sim.Nr);
BFlux = zeros(length(sim.t),30);


for i = 1:length(sim.t)

    % dry mass per volume in each time step in each size-density bin [µgC/m3]
    Mdry(i,:,:) = reshape(squeeze(sim.Mo(i,:)),sim.Nd,sim.Nr);


    % flux out of the H=50m mixed layer --> divide /H ???  [µgC/m2/day]
    Flux(i,:,:) = squeeze(Mdry(i,:,:)).*sim.w;

    % Size integrated flux [µgC/m2/day]
    BFlux(i,:) = sum(squeeze(Flux(i,:,:)),1);

end

%%

% real size values:
% r = sim.p.ro*sim.p.delta.^sim.x;
% dr = r*(1-1/sim.p.delta);

figure;
contourf(sim.r_mean, sim.t, BFlux)
set(gca, 'XScale', 'log')
xlabel('Radius (µm)')
ylabel('Time [days]')
c = colorbar;
c.Label.String = 'Size-integr. Flux [µgC/m2/day]';
c.Label.FontSize = 11;
% title = 'Size-integr. Flux [µgC/m2/day]';
% xlim([min(r) 8.2536e+04])
% xlim([1 8.2536e+04])
xlim([0.1946 7.9552e+05])

%% Productivity

pr = zeros(length(sim.t),sim.Nd,sim.Nr);
pr2 = zeros(length(sim.t),30);
tt = sim.t;

for i = 1:length(sim.t)

    pr(i,:,:) = sim.prod*(1-cos(2*pi*i/365))/2;

    pr2(i,:) = sum(squeeze(pr(i,:,:)),1);


end

%%

figure;
contourf(sim.x(1:5), sim.t, pr2(:,1:5))
xlabel('X')
ylabel('Time [days]')
c = colorbar;
c.Label.String = 'Productivity [µgC/m2/day]';
c.Label.FontSize = 11;

%%

figure;
pr3 = sum(pr(end-365:end,:,1),2);
plot(sim.t(end-365:end), pr3);


%%
% 4800 gC/ m2/ day
% 72 ^6
% 26 ^5

vv = sum(sum(BFlux(end-365:end,:),2)) / 10^6
%%

BFlux_1year = sum(sum(BFlux(end-365:end,:)))*1E-6 % [gC / m2/ year] 72 vs  26


