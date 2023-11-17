option = 9;

if option == 1
    % param = d_mean*1E-18*1E9;
    param = d_mean;
    r_sec = r_mean;
    % t_text = 'Excess density \rho [\mug / \mum^3]';
    t_text = 'Excess density \rho [kg / m^3]';
elseif option == 2
    param = phi;
    r_sec = r_mean;
    t_text = 'Porosity estimate [-]';

elseif  option == 3
    param = m_mean(:,1:16);
    r_sec = r_mean(1:16);
    t_text = 'Total mass [\mug]';

elseif option == 4
    param = mdry2(:,1:10);
    r_sec = r_mean(1:10);
    t_text = 'Dry mass NEW [\mug]';


elseif option == 5
    param = mdry(:,1:10);
    r_sec = r_mean(1:10);
    t_text = 'Dry mass OLD [\mug]';

elseif option == 6
    param = mdry(:,1:16)-mdry2(:,1:16);
    r_sec = r_mean(1:16);
    t_text = 'Dry mass DIFF (old-new) [\mug]';

elseif option == 7
    param = w(:,1:21);
    r_sec = r_mean(1:21);
    t_text = 'Sinking velocity [m/day]';

elseif option == 8
    param = pfrag;
    r_sec = r_mean;
    t_text = 'Fragmentation rate [1/day]';

elseif option == 9
    % remin = Rrate.*remin_grid;
    param = remin;
    r_sec = r_mean;
    t_text = 'Remineralization rate [1/day]';

elseif option == 10
    param = m_mean;
    r_sec = r_mean;
    t_text = 'Total mass [\mug]';

elseif option == 11
    param = reminDif;
    % r_sec = r_mean;
    t_text = 'Remineralization rate [1/day]';

end

% d_real = drho*z * 1E-9;
d_real = rhoo*lambda.^z;

figure;
imagesc(r_sec, d_real, param);
% imagesc(param);
colorbar
set(gca,'ColorScale','log', 'XScale', 'log', 'YScale', 'log')
% set(gca,'ColorScale','log')
% clim([0.1 200])

axis xy
colorbar
% set(gca,'ColorScale','log', 'XScale', 'log')
% set(gca,'XScale', 'log')

title(t_text)
xlabel('radius [\mum]')
ylabel('density interval [kg / m^3]')


%% delta

dd = q.^x;

semilogx(r_mean, dd, 'LineWidth', 2)
xlabel('radius [\mum]')
ylabel('porosity estimate (-)')
xlim([min(r_mean) max(r_mean)])