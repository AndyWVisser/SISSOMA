% M = Mo(end,:);
% m = mdry;

dNrem2 = remin(:).*M(:);
dNrom2 = -dNrem2 + [dNrem2(2:end); 0]; % the last part corresponds to a shift of 1 to have the Z+1 * N(Z+1) 
dNrom2(Nd:Nd:end) = -dNrem2(Nd:Nd:end);
dMremin2 = dNrom2(:);  % .*m(:); % go back to dM (mass multiplying by mass per particle)
dMremin2(1:Nd:end) = dMremin2(1:Nd:end)- Rrate*M(1:Nd:end);  % why add the Rrate*M ??? --> ensure remin make mass leave syst at Z=1



dNrem = remin(:).*N(:);
dNrom = -dNrem + [dNrem(2:end); 0]; % the last part corresponds to a shift of 1 to have the Z+1 * N(Z+1) 
dNrom(Nd:Nd:end) = -dNrem(Nd:Nd:end);
dMremin = dNrom(:).*m(:); % go back to dM (mass multiplying by mass per particle)
dMremin(1:Nd:end) = dMremin(1:Nd:end)- Rrate*M(1:Nd:end);  % why add the Rrate*M ??? --> ensure remin make mass leave syst at Z=1


reminDif = dMremin - dMremin2;


% dNrem = remin(:).*N(:);
% dNrom = -dNrem + [dNrem(2:end); 0]; % the last part corresponds to a shift of 1 to have the Z+1 * N(Z+1) 
% dNrom(Nd:Nd:end) = -dNrem(Nd:Nd:end);
% dMremin = dNrom(:).*m(:); % go back to dM (mass multiplying by mass per particle)
% dMremin(1:Nd:end) = dMremin(1:Nd:end)- Rrate*M(1:Nd:end);  % why add the Rrate*M ??? --> ensure remin make mass leave syst at Z=1