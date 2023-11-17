% figure;

% semilogx(sim.r_mean,sim.BFlux,'-o','linewidth', 1.1); 
plot(sim.x+0.5,sim.BFlux,'-o','linewidth', 1.1); 


title('size integrated flux [µgC/m2/day]');
hold on 
% xlim([0.25 1E6])
% xlabel('Size (µm)')
xlabel('X')

ylabel('size integrated flux [µgC/m2/day]')

%legend('min=1, max=1E6', 'min=1, max=1E5', 'min=0.5, max=1E6' , 'min=0.5, max=1E5', 'min=0.25, max=1E6', 'min=0.25, max=1E5', 'min=0.25, max=1E4' )


