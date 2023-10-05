simA = load("sim_scenario_A.mat");
simA = simA.sim;

maxFLUX_A = max(max(simA.Flux));
maxDrymass_A = max(max(simA.Mdry));
maxBFLUX_A = max(max(simA.BFlux));
maxdMdt_A = max(max(simA.dMdto));
maxdMsink_A = max(max(simA.dMsink));
maxdMfrag_A = max(max(simA.dMfrag));
maxdMremin_A = max(max(simA.dMremin));
maxdMaggre_A = max(max(simA.dMaggreg));


minFLUX_A = min(min(simA.Flux));
minDrymass_A = min(min(simA.Mdry));
minBFLUX_A = min(min(simA.BFlux));
mindMdt_A = min(min(simA.dMdto));
mindMsink_A = min(min(simA.dMsink));
mindMfrag_A = min(min(simA.dMfrag));
mindMremin_A = min(min(simA.dMremin));
mindMaggre_A = min(min(simA.dMaggreg));

%%

simB = load("sim_scenario_B.mat");
simB = simB.sim;

maxFLUX_B = max(max(simB.Flux));
maxDrymass_B = max(max(simB.Mdry));
maxBFLUX_B = max(max(simB.BFlux));
maxdMdt_B = max(max(simB.dMdto));
maxdMsink_B = max(max(simB.dMsink));
maxdMfrag_B = max(max(simB.dMfrag));
maxdMremin_B = max(max(simB.dMremin));
maxdMaggre_B = max(max(simB.dMaggreg));


minFLUX_B = min(min(simB.Flux));
minDrymass_B = min(min(simB.Mdry));
minBFLUX_B = min(min(simB.BFlux));
mindMdt_B = min(min(simB.dMdto));
mindMsink_B = min(min(simB.dMsink));
mindMfrag_B = min(min(simB.dMfrag));
mindMremin_B = min(min(simB.dMremin));
mindMaggre_B = min(min(simB.dMaggreg));

%% 

maxFLUX = max(maxFLUX_A, maxFLUX_B);
maxDrymass = max(maxDrymass_A, maxDrymass_B);
maxBFLUX = max(maxBFLUX_A, maxBFLUX_B);
maxdMdt = max(maxdMdt_A, maxdMdt_B);
maxdMsink = max(maxdMsink_A, maxdMsink_B);
maxdMfrag = max(maxdMfrag_A, maxdMfrag_B);
maxdMremin = max(maxdMremin_A, maxdMremin_B);
maxdMaggre = max(maxdMaggre_A, maxdMaggre_B);


minFLUX = min(minFLUX_A, minFLUX_B);
minDrymass = min(minDrymass_A, minDrymass_B);
minBFLUX = min(minBFLUX_A, minBFLUX_B);
mindMdt = min(mindMdt_A, mindMdt_B);
mindMsink = min(mindMsink_A, mindMsink_B);
mindMfrag = min(mindMfrag_A, mindMfrag_B);
mindMremin = min(mindMremin_A, mindMremin_B);
mindMaggre = min(mindMaggre_A, mindMaggre_B);


%% PACKAGE 1

%scenario A:
x = simA.x;
z = simA.z;
Mdry = simA.Mdry;
Flux = simA.Flux;
BFlux = simA.BFlux;
RMSE = simA.RMSE;
t = simA.t;

figure; 
subplot(2,2,1); imagesc(x,z,Mdry); colorbar; clim([minDrymass maxDrymass]) ; title('dry mass'); axis xy
subplot(2,2,2); imagesc(x,z,Flux); colorbar; clim([minFLUX maxFLUX]); title('flux'); axis xy
subplot(2,2,4); semilogy(t(2:end),RMSE,'.');  title('root mean error');
subplot(2,2,3); plot(x+.5,BFlux,'o'); ylim([minBFLUX maxBFLUX]); title('size integrated flux');

%%

%scenario B:
x = simB.x;
z = simB.z;
Mdry = simB.Mdry;
Flux = simB.Flux;
BFlux = simB.BFlux;
RMSE = simB.RMSE;
t = simB.t;

figure; 
subplot(2,2,1); imagesc(x,z,Mdry); colorbar; clim([minDrymass maxDrymass]) ; title('dry mass'); axis xy
subplot(2,2,2); imagesc(x,z,Flux); colorbar; clim([minFLUX maxFLUX]); title('flux'); axis xy
subplot(2,2,4); semilogy(t(2:end),RMSE,'.');  title('root mean error');
subplot(2,2,3); plot(x+.5,BFlux,'o'); ylim([minBFLUX maxBFLUX]); title('size integrated flux');

%% diff:

figure; 
subplot(2,2,1); imagesc(x,z,simA.Mdry - simB.Mdry); colorbar; title('dry mass'); axis xy
subplot(2,2,2); imagesc(x,z,simA.Flux - simB.Flux); colorbar; title('flux'); axis xy
% subplot(2,2,4); semilogy(t(2:end),RMSE,'.');  title('root mean error');
subplot(2,2,3:4); plot(x+.5,simA.BFlux,'o', 'Linewidth', 1.5); hold on ; plot(x+.5,simB.BFlux,'x', 'Linewidth', 1.5); legend('scenario A', 'scenatrio B') ; title('size integrated flux');



%% PACKAGE 2

% scenario A:

dMdto = simA.dMdto;
dMsink = simA.dMsink;
dMremin = simA.dMremin;
dMfrag = simA.dMfrag;
dMaggreg = simA.dMaggreg;
prod = simA.prod;

Nd = simA.Nd;
Nr = simA.Nr;



figure; 
subplot(3,2,1); imagesc(x,z,reshape(dMdto,Nd,Nr)); colorbar; clim([mindMdt maxdMdt]) ; title('dMdt'); axis xy;
subplot(3,2,2); imagesc(x,z,reshape(dMsink,Nd,Nr)); colorbar; clim([mindMsink maxdMsink]) ; title('sinking'); axis xy;
subplot(3,2,4); imagesc(x,z,reshape(dMremin,Nd,Nr)); colorbar; clim([mindMremin maxdMremin]) ; title('remineralization'); axis xy;
subplot(3,2,3); imagesc(x,z,reshape(dMfrag,Nd,Nr)); colorbar; clim([mindMfrag maxdMfrag]) ; title('fragmentation'); axis xy;
subplot(3,2,5); imagesc(x,z,prod); colorbar; title('production'); axis xy;
subplot(3,2,6); imagesc(x,z,reshape(dMaggreg,Nd,Nr)); colorbar; clim([mindMaggre maxdMaggre]) ; title('aggregation'); axis xy;

%%
% scenario B:

dMdto = simB.dMdto;
dMsink = simB.dMsink;
dMremin = simB.dMremin;
dMfrag = simB.dMfrag;
dMaggreg = simB.dMaggreg;
prod = simB.prod;

Nd = simA.Nd;
Nr = simA.Nr;



figure; 
subplot(3,2,1); imagesc(x,z,reshape(dMdto,Nd,Nr)); colorbar; clim([mindMdt maxdMdt]) ; title('dMdt'); axis xy;
subplot(3,2,2); imagesc(x,z,reshape(dMsink,Nd,Nr)); colorbar; clim([mindMsink maxdMsink]) ; title('sinking'); axis xy;
subplot(3,2,4); imagesc(x,z,reshape(dMremin,Nd,Nr)); colorbar; clim([mindMremin maxdMremin]) ; title('remineralization'); axis xy;
subplot(3,2,3); imagesc(x,z,reshape(dMfrag,Nd,Nr)); colorbar; clim([mindMfrag maxdMfrag]) ; title('fragmentation'); axis xy;
subplot(3,2,5); imagesc(x,z,prod); colorbar; title('production'); axis xy;
subplot(3,2,6); imagesc(x,z,reshape(dMaggreg,Nd,Nr)); colorbar; clim([mindMaggre maxdMaggre]) ; title('aggregation'); axis xy;

%% diff:

figure; 
subplot(3,2,1); imagesc(x,z,reshape(simA.dMdto - simB.dMdto,Nd,Nr)); colorbar; title('dMdt'); axis xy;
subplot(3,2,2); imagesc(x,z,reshape(simA.dMsink - simB.dMsink,Nd,Nr)); colorbar; title('sinking'); axis xy;
subplot(3,2,4); imagesc(x,z,reshape(simA.dMremin - simB.dMremin,Nd,Nr)); colorbar; title('remineralization'); axis xy;
subplot(3,2,3); imagesc(x,z,reshape(simA.dMfrag - simB.dMfrag,Nd,Nr)); colorbar; title('fragmentation'); axis xy;
subplot(3,2,5); imagesc(x,z,prod); colorbar; title('production'); axis xy;
subplot(3,2,6); imagesc(x,z,reshape(simA.dMaggreg - simB.dMaggreg,Nd,Nr)); colorbar; title('aggregation'); axis xy;










