Temps = 7:13;

a = 2;           % selfsimilarity parameter; typically between 1.8 and 2.0
alpha = 0.4;      % stickiness; range 0 to 1
epsilon = 1E-6;  % turbulent dissipation rate [m3/s2]
Ptotal = 1E6;    % total productivity; typically 1E6 [µg m-2 day-1] (1 gC m-2 day-1)
Rrate = 0.1;     % remineralization rate [day-1]
Frate = 500;     % maximum fragmentation rate [day-1] for aggregates > 1 m
Tmax = 5*365;    % period of simulation [days]
seasonal = 0;   % seasonal (true) or constant (false) production
temp_depend_remin = 1; % constant remineralization rate (false) or temperature dependent (true)
T_input = 10;


Bfluxs = zeros(length(Temps), 30);
Rates = zeros(length(Temps), 1);



for i = 1:length(Temps)

    sim = coagfunDev(a,alpha,epsilon,Ptotal,Rrate,Frate,Tmax,seasonal,temp_depend_remin,Temps(i))

    % dry mass per volume in last time step in each size-density bin [µgC/m3]
    Mdry = reshape(sim.Mo(end,:),sim.Nd,sim.Nr);

    % flux out of the H=50m mixed layer --> divide /H ???     [µgC/m2/day]
    Flux = Mdry.*sim.w;

    %
    BFlux = sum(Flux,1);


    Bfluxs(i,:) = BFlux;
    Rates(i,1) = sim.Rrate;

end


%%

figure

for j = 1:length(Temps)

plot(sim.x+.5, Bfluxs(j,:)','-', 'LineWidth', 1.5); 

hold on
end

xlabel('r [\mum]')
ylabel('size integrated flux [µgC/m2/day]');
leg = legend('7', '8', '9', '10', '11', '12', '13')
title(leg, 'Temperature (C)')


%%
figure

for j = 1:length(Temps)

plot(Rates(j,1), 'o', 'LineWidth', 2)

hold on
end

% xlabel('r [\mum]')
ylabel('Rrate [day-1]');
leg = legend('7', '8', '9', '10', '11', '12', '13');
title(leg, 'Temperature (C)')
