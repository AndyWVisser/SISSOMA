function sensitivity_analysis3_seasonal_f_normP(param_number)

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
Tmax  = 5*365;     % period of simulation [days]
seasonal = 1;

myCols = {[0.2 0.5 0.9 1]   [0.7 0.5 0.9 1] [0.2 0.1 0.2 1]  [0.6 0.9 0.3 1]};
myShades = {[0.2 0.5 0.9 0.08] [0.7 0.5 0.9 0.08] [0.2 0.1 0.2 0.08]  [0.6 0.9 0.3 0.08]};

figure;

if param_number == 1
    %     a_new = [1.7 1.8 1.9 2.0 2.1 2.2 2.3];
    a_new = [1.7 1.9 2.1 2.3];
    legendtext = {'a = '};
    legendtext0 = {'1.7','1.9','2.1','2.3'};
    curr_analy = a_new;
    code = 'Self-similarity parameter [-]';
    code_leg = '\alpha (self-similarity) ';

elseif param_number == 2
    alpha_new = [0.1 0.4 0.7 1]; % [0.3, 0.5, 0.7, 0.9];
    legendtext = {'alpha = '};
    legendtext0 = {'0.1','0.4','0.7','1'}; % {'0.3','0.5','0.7','0.9'};
    curr_analy = alpha_new;
    code = 'Stickiness [-]';
    code_leg = 'a (stickiness)';

elseif param_number == 3
    epsilon_new = [1E-2, 1E-4, 1E-6 1E-8];
    legendtext = {'epsilon = '};
    legendtext0 = {'10^{-2}','10^{-4}','10^{-6}','10^{-8}'};
    curr_analy = epsilon_new;
    code = 'Turbulent dissipation rate [m^{3}/s^{2}]';
    code_leg = '\epsilon';

elseif param_number == 4
    P_new = [1E3 1E4 1E5 1E6];
    legendtext = {'Ptot = '};
    legendtext0 = {'10^{3}', '10^{4}','10^{5}','10^{6}'};
    curr_analy = P_new;
    code = 'Total productivity [{\mu}gC/ m^{2}/ day]';
    code_leg = 'P';

end

close all
for i = 1:length(curr_analy)

    % Frate = Frate*Ptotal*1E-6;

    switch param_number
        case 1
            sim = coagfunDev(a_new(i),alpha,epsilon,Ptotal,Rrate,Frate,Tmax,seasonal);

        case 2
            sim = coagfunDev(a,alpha_new(i),epsilon,Ptotal,Rrate,Frate,Tmax,seasonal);

        case 3
            sim = coagfunDev(a,alpha,epsilon_new(i),Ptotal,Rrate,Frate,Tmax,seasonal);

        case 4
            %             Frate = Frate*P_new(i)*1E-6;
            sim = coagfunDev(a,alpha,epsilon,P_new(i),Rrate,Frate,Tmax,seasonal);
    end


    LY = length(sim.t)-366:length(sim.t);

    Mdry = zeros(length(LY), sim.Nd, sim.Nr);
    Flux = zeros(length(LY), sim.Nd, sim.Nr);
    BFlux = zeros(length(LY), sim.Nr);

    % Get last year
    for iTime = 1:length(LY)

        % dry mass per volume in each time step in each size-density bin [µgC/m3]
        Mdry(iTime,:,:) = reshape(squeeze(sim.Mo(LY(iTime),:)),sim.Nd,sim.Nr);


        % flux out of the H=50m mixed layer  [µgC/m2/day]
        Flux(iTime,:,:) = squeeze(Mdry(iTime,:,:)).*sim.w;

        fseas = (1-cos(2*pi*LY(iTime)/365))/2.5 + 0.2;
        
        % RATIO: Size integrated flux [µgC/m2/day] / Production in each time step
        BFlux(iTime,:) = sum(squeeze(Flux(iTime,:,:)),1) ./ (sim.p.Ptotal*fseas);

    end

    for k = 1:length(LY)

        semilogx(sim.r_mean, BFlux(k,:), 'LineWidth', 0.2, 'Color', myShades{i})
        % loglog(sim.r_mean, BFlux(k,:), 'LineWidth', 0.2, 'Color', myShades{i})
        hold on
    end

    Bmean = mean(BFlux,1);
    myMean(i) = semilogx(sim.r_mean, Bmean, '-o', 'LineWidth', 1, 'Color', myCols{i});
    % myMean(i) = loglog(sim.r_mean, Bmean, '-o', 'LineWidth', 1, 'Color', myCols{i});

end

handlevec = [myMean(1) myMean(2) myMean(3) myMean(4)];
leg = legend(handlevec, legendtext0);
title(leg,code_leg)

xlabel('r [\mum]')
ylabel('{\it f}(r)', 'FontSize', 17)

xlim([min(sim.r_mean) max(sim.r_mean)])

% if param_number == 4
    % ylim([10^-2 inf])
% end
end





