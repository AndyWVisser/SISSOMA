function basic_plots(sim)

x = sim.x;
z = sim.z;
Mdry = sim.Mdry;
Flux = sim.Flux;
BFlux = sim.BFlux;
RMSE = sim.RMSE;
t = sim.t;

% prod = sim.prod;
% dMdto = sim.dMdto;
% dMsink = sim.dMsink;
% dMremin = sim.dMremin;
% dMfrag = sim.dMfrag;
% dMaggreg = sim.dMaggreg;
% Nd = sim.Nd;
% Nr = sim.Nr;

%%
figure(1); clf;
subplot(2,2,1); imagesc(x,z,Mdry); colorbar;  title('dry mass'); axis xy
subplot(2,2,2); imagesc(x,z,Flux); colorbar; title('flux'); axis xy
subplot(2,2,4); semilogy(t(2:end),RMSE,'.'); title('root mean error');
subplot(2,2,3); plot(x+.5,BFlux,'o'); title('size integrated flux');

%%
% figure(2); clf;
% subplot(3,2,1); imagesc(x,z,reshape(dMdto,Nd,Nr)); colorbar; title('dMdt'); axis xy;
% subplot(3,2,2); imagesc(x,z,reshape(dMsink,Nd,Nr)); colorbar; title('sinking'); axis xy;
% subplot(3,2,4); imagesc(x,z,reshape(dMremin,Nd,Nr)); colorbar; title('remineralization'); axis xy;
% subplot(3,2,3); imagesc(x,z,reshape(dMfrag,Nd,Nr)); colorbar; title('fragmentation'); axis xy;
% subplot(3,2,5); imagesc(x,z,prod); colorbar; title('production'); axis xy;
% subplot(3,2,6); imagesc(x,z,reshape(dMaggreg,Nd,Nr)); colorbar; title('aggregation'); axis xy;


end