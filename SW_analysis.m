% T_input = 2:4:22;
lowlim = 15;
uplim = 20;

T_input = [lowlim lowlim linspace(lowlim,uplim,182) linspace(uplim,lowlim,182) linspace(lowlim,uplim,182) linspace(uplim,lowlim,182) linspace(lowlim,uplim,182) linspace(uplim,lowlim,182) linspace(lowlim,uplim,182) linspace(uplim,lowlim,182) linspace(lowlim,uplim,182) linspace(uplim,lowlim,182) lowlim lowlim lowlim ];


rho_sw = zeros(length(T_input),1);
nu = zeros(length(T_input),1);
mu = zeros(length(T_input),1);

for j = 1:length(T_input)
    rho_sw(j) =  SW_Density(T_input(j), 'C', 35, 'ppt', 3, 'bar');
    nu(j) = SW_Kviscosity(T_input(j), 'C', 35, 'ppt');
    mu(j) = nu(j)*rho_sw(j);

end

rho_sw_M = mean(rho_sw);
nu_M = mean(nu);
mu_M = mean(mu);
%%

col = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330] };

x = 1:1825;
% figure;
plot(x,mu, '-', 'LineWidth', 1, 'Color', col{4})
hold on
plot(x, ones(1,1825)*mu_M, '--', 'LineWidth', 1.2, 'Color', col{4})
xlabel('Time [days]')
% ylabel('Water density [kg m^-3]')
% ylabel('Kinematic viscosity [[m^2/ s]]')
ylabel('Dynamic viscosity [kg/ m/ s]')
leg = legend('0-5', '5-10', '10-15', '15-20');
title(leg, 'Temp range')

%% Salinity!!!!!!

salins = linspace(33,37, 10);
temps = linspace(5,20,1);

figure;

for iS = 1:length(salins)
    for iT = 1:length(temps)

        rho_sw =  SW_Density(temps(iT), 'C', salins(iS), 'ppt', 3, 'bar');
        nu = SW_Kviscosity(temps(iT), 'C', salins(iS), 'ppt');
        mu = nu*rho_sw;

        hold on
        % plot(rho_sw, 'o')
        % plot(nu, 'o')
        plot(mu, 'o')

    end
end


xlabel('Time [days]')
% ylabel('Water density [kg m^-3]')
% ylabel('Kinematic viscosity [m^2/ s]')
ylabel('Dynamic viscosity [kg/ m/ s]')










