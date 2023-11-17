figure
semilogx(r_mean, mean(w), 'LineWidth', 1.5)
hold on
semilogx(r_mean, max(w), '--', 'LineWidth', 1.4)
semilogx(r_mean, min(w), '--', 'LineWidth', 1.4)

legend('mean', 'max', 'min')

xlabel('X')
ylabel('Z')
% colorbar
title('sinking speed [m/ day]')

%% 
figure
imagesc(x, z, Re); axis xy
colorbar