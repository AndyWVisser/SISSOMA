clear all

%coagfun001(a,alpha,epsilon,Ptot,Rrate,Frate,Tmax,Season);

a = 2.0;
alpha = 0.2;
epsilon = 1E-6
Ptot = 1E6
Rrate = 0.1;
Frate = 100;
Tmax = 2000;
seasonal = 0;

simf010 = coagfun001(a,alpha,epsilon,Ptot,Rrate,0.1*Frate,Tmax,seasonal);
simf050 = coagfun001(a,alpha,epsilon,Ptot,Rrate,0.5*Frate,Tmax,seasonal);
simf100 = coagfun001(a,alpha,epsilon,Ptot,Rrate,Frate,Tmax,seasonal);
simf200 = coagfun001(a,alpha,epsilon,Ptot,Rrate,2*Frate,Tmax,seasonal);
simf500 = coagfun001(a,alpha,epsilon,Ptot,Rrate,5*Frate,Tmax,seasonal);

% sima180 = coagfun001(1.80,alpha,epsilon,Ptot,Rrate,Frate,Tmax,seasonal);
% sima185 = coagfun001(1.85,alpha,epsilon,Ptot,Rrate,Frate,Tmax,seasonal);
sima190 = coagfun001(1.90,alpha,epsilon,Ptot,Rrate,Frate,Tmax,seasonal);
sima195 = coagfun001(1.95,alpha,epsilon,Ptot,Rrate,Frate,Tmax,seasonal);
sima200 = coagfun001(a,alpha,epsilon,Ptot,Rrate,Frate,Tmax,seasonal);
sima205 = coagfun001(2.05,alpha,epsilon,Ptot,Rrate,Frate,Tmax,seasonal);
sima210 = coagfun001(2.1,alpha,epsilon,Ptot,Rrate,Frate,Tmax,seasonal);

sims01 = coagfun001(a,0.5*alpha,epsilon,Ptot,Rrate,Frate,Tmax,seasonal);
sims02 = coagfun001(a,0.75*alpha,epsilon,Ptot,Rrate,Frate,Tmax,seasonal);
sims03 = coagfun001(a,alpha,epsilon,Ptot,Rrate,Frate,Tmax,seasonal);
sims04 = coagfun001(a,1.25*alpha,epsilon,Ptot,Rrate,Frate,Tmax,seasonal);
sims05 = coagfun001(a,1.5*alpha,epsilon,Ptot,Rrate,Frate,Tmax,seasonal);
%sims06 = coagfun001(a,2*alpha,epsilon,Ptot,Rrate,Frate,Tmax,seasonal);


sime08 = coagfun001(a,alpha,0.01*epsilon,Ptot,Rrate,Frate,Tmax,seasonal);
sime07 = coagfun001(a,alpha,0.1*epsilon,Ptot,Rrate,Frate,Tmax,seasonal);
sime06 = coagfun001(a,alpha,epsilon,Ptot,Rrate,Frate,Tmax,seasonal);
sime05 = coagfun001(a,alpha,10*epsilon,Ptot,Rrate,Frate,Tmax,seasonal);
sime04 = coagfun001(a,alpha,100*epsilon,Ptot,Rrate,Frate,Tmax,seasonal);
%sime03 = coagfun001(a,alpha,1000*epsilon,Ptot,Rrate,Frate,Tmax,seasonal);

%simp01 = coagfun001(a,alpha,epsilon,0.1*Ptot,Rrate,Frate,Tmax,seasonal);
simp02 = coagfun001(a,alpha,epsilon,0.2*Ptot,Rrate,Frate,Tmax,seasonal);
simp05 = coagfun001(a,alpha,epsilon,0.5*Ptot,Rrate,Frate,Tmax,seasonal);
simp10 = coagfun001(a,alpha,epsilon,Ptot,Rrate,Frate,Tmax,seasonal);
simp20 = coagfun001(a,alpha,epsilon,2*Ptot,Rrate,Frate,Tmax,seasonal);
simp50 = coagfun001(a,alpha,epsilon,5*Ptot,Rrate,Frate,Tmax,seasonal);
%simp100 = coagfun001(a,alpha,epsilon,10*Ptot,Rrate,Frate,Tmax,seasonal);

save simsall001;
