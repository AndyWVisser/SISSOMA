%%

dNrem = remin(:).*N(:);
dNrom = -dNrem + [dNrem(2:end); 0]; % the last part corresponds to a shift of 1 to have the Z+1 * N(Z+1) 
dNrom(Nd:Nd:end) = -dNrem(Nd:Nd:end); % fix for bins: 10,20,30... (high density)--> nothing comes in
dMremin = dNrom(:).*m(:); % go back to dM (mass multiplying by mass per particle)
dMremin(1:Nd:end) = dMremin(1:Nd:end) - Rrate*M(1:Nd:end);  % ensure remin make mass to leave the system at Z=1 (low density)


dMrem2 = remin(:) .* M(:);
dMremin2 = -dMrem2 + [dMrem2(2:end); 0];
dMremin2(Nd:Nd:end) = -dMrem2(Nd:Nd:end);
dMremin2(1:Nd:end) = dMremin2(1:Nd:end) - Rrate*M(1:Nd:end);



%%
CC = (dMremin - dMremin2);
c = CC./dMremin;


figure
imagesc(Nr,Nd,reshape(CC,Nd,Nr))
axis xy
colorbar


figure
imagesc(Nr,Nd,reshape(CC./dMremin,Nd,Nr))
axis xy
colorbar

figure
rat = dMremin./dMremin2;
rat(isnan(rat)) = 1
imagesc(Nr,Nd,reshape(rat,Nd,Nr))
axis xy
colorbar
title('scenario_N / scenario_M')

