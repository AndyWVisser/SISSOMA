figure

simP = load('sim_PARTICLES.mat');
simP = simP.sim;
simM = load('sim_MASS.mat');
simM = simM.sim;

x = simP.x;
z = simP.z;
t = simP.t;


simP.Mdry(simP.Mdry < 0.01) = 0.01;
simM.Mdry(simM.Mdry < 0.01) = 0.01;

simP.Flux(simP.Flux < 0.01) = 0.01;
simM.Flux(simM.Flux < 0.01) = 0.01;


subplot(2,2,1); imagesc(x,z,simP.Mdry./ simM.Mdry); colorbar;  title('dry mass'); axis xy ; % clim([0 1]);
subplot(2,2,2); imagesc(x,z,simP.Flux./simM.Flux); colorbar; title('flux'); axis xy
% subplot(2,2,4); semilogy(t(2:end),RMSE,'.'); title('root mean error');
subplot(2,2,3:4); plot(x+.5,simP.BFlux,'o'); title('size integrated flux'); hold on ; plot(x+.5,simM.BFlux,'x') ; legend('scenario-Particles', 'scenario-Mass')