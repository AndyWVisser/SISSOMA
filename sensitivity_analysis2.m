function sensitivity_analysis2(param_number)

% param_number = 1:
a = 2;            % selfsimilarity parameter; typically between 1.8 and 2.0
% param_number = 2:
alpha = 0.3;      % stickiness; range 0 to 1
% param_number = 3:
epsilon = 1E-7;   % turbulent dissipation rate [m2/s3]
% param_number = 4:
Ptotal = 1E6;     % total productivity; typically 1E6 [µg  m-2 day-1] (1 gC m-2 day-1)

Frate = 500;      % maximum fragmentation rate [day-1] for aggregates > 1 m
Rrate  = 0.1;     % remineralization rate [day-1]
Tmax  = 2000;     % period of simulation [days]
seasonal = 0;

if param_number == 1
    a_new = [1.7 1.8 1.9 2.0 2.1 2.2 2.3];
    legendtext = {'a = '};
    legendtext0 = {'1.7','1.8','1.9','2.0','2.1','2.2','2.3'};
    curr_analy = a_new;
    code = 'Sef-similarity parameter [-]';
    code_leg = '\alpha';

elseif param_number == 2
    alpha_new = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7];
    legendtext = {'alpha = '};
    legendtext0 = {'0.1','0.2','0.3','0.4','0.5','0.6','0.7'};
    curr_analy = alpha_new;
    code = 'Stickiness [-]';
    code_leg = 'a';

elseif param_number == 3
    epsilon_new = [1E-2, 1E-3, 1E-4, 1E-5, 1E-6, 1E-7, 1E-8];
    legendtext = {'epsilon = '};
    legendtext0 = {'10^{-2}','10^{-3}','10^{-4}','10^{-5}','10^{-6}','10^{-7}','10^{-8}'};
    curr_analy = epsilon_new;
    code = 'Turbulent dissipation rate [m^{3}/s^{2}]';
    code_leg = '\epsilon';

elseif param_number == 4
    P_new = [1E5 1E6 1E7 1E8 1E9];
    legendtext = {'Ptot = '};
    legendtext0 = {'10^{5}','10^{6}','10^{7}','10^{8}','10^{9}'};
    curr_analy = P_new;
    code = 'Total productivity [{\mu}gC/ m^{2}/ day]';
    code_leg = 'P';

end

close all
for i = 1:length(curr_analy)

    Frate = Frate*Ptotal*1E-6;
   
    switch param_number
        case 1
            sim = coagfunDev(a_new(i),alpha,epsilon,Ptotal,Rrate,Frate,Tmax,seasonal);

        case 2
            sim = coagfunDev(a,alpha_new(i),epsilon,Ptotal,Rrate,Frate,Tmax,seasonal);

        case 3
            sim = coagfunDev(a,alpha,epsilon_new(i),Ptotal,Rrate,Frate,Tmax,seasonal);

        case 4
            Frate = Frate*P_new(i)*1E-6;
            sim = coagfun002(a,alpha,epsilon,P_new(i),Rrate,Frate,Tmax,seasonal);
    end

    %% Plot
    maxValue = -inf;

    % Upper limit??:
    r = sim.p.ro*sim.p.delta.^sim.x;

    % calculate dr:
    dr = r*(1-1/sim.p.delta);

    % Normalized number of particles/ size bin?:
    N_norm = sum(sim.N)./dr;
    % N_norm = cumsum(sum(sim.N)./dr,"reverse");

    % Normalized Flux (/Ptotal):
    if param_number == 4
        Flux_norm = sim.BFlux / curr_analy(i);
    else
        Flux_norm = sim.BFlux / Ptotal;
    end

    m = sim.x + 0.5;
    curr_value = curr_analy(i)*ones(1,30);


    figure(1)
    semilogx(r,Flux_norm,'LineWidth',1.3);
    ag = gca; ag.FontSize = 11;
    xlabel('r [\mum]')
    ylabel('{\it f}(r)', 'FontSize', 17)
    % title(code)
    hold on

    figure(2)
    loglog(r,N_norm,'LineWidth',1.3);
    an = gca; an.FontSize = 11;
    ylim([1E-10 1E15])
    xlim([1E0 1E6])
    % title(code)
    xlabel('r [\mum]')
    ylabel('{\it n}(r)', 'FontSize', 17)
    hold on


end
% legend_fin = strcat(legendtext,legendtext0);
% legend(ag, legendtext0)

figure(2)
hold on

s1 = 1E13*r.^(-3);
s2 = 1E13*r.^(-4);

loglog(r,s1,'--','Color', [0.7,0.7,0.7],'LineWidth',1);
loglog(r,s2,'--', 'Color', [0,0,0], 'LineWidth',1);

legendtext0{end+1} = 'Slope of -3';
legendtext0{end+1} = 'Slope of -4';
leg = legend(legendtext0);
title(leg, code_leg)
end






