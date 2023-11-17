col = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330] };

figure;

loglog(r_mean, max(d_mean), 'LineWidth', 1.2, 'Color', col{1})
xlim([0.25 1E6])
hold on
loglog(r_mean, min(d_mean), '--', 'LineWidth', 1.2, 'Color', col{1})

xlabel('radius [\mum]')
ylabel('Excess density[kg m-3]')
title('Excess density[kg m-3]')
% legend('min=1, max=1E6', 'min=0.25, max=1E5')
