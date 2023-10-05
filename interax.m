function [dMdt,dMsink,dMremin,dMfrag,dMaggr] = interax(t,M,m,bi,bj,Nr,Nd,b300,b301,b310,b311,f00,f01,f10,f11,alpha,beta,w,H,prod,remin,Rrate,pfrag)
% function dMdt = interax(t,M,m,bi,bj,Nr,Nd,b300,b301,b310,b311,f00,f01,f10,f11,alpha,beta,w,H,prod,remin,Rrate,pfrag)

% Do we calculate dry mass or total(dry+water)?

% Dry mass?? [µgC/m3]
M(M<0) = 0;

% Dry mass per particle reference
m = m(:); 

% Number of particles per volume [#/ m3]
N = M(:)./m(:);

dMaggr = zeros(size(M));
% dMremin = zeros(size(dMaggr));
% dMfrag = zeros(size(dMaggr));


%% Aggregation
for k = 1:length(bi)
    ii = bi(k)+1;
    jj = bj(k)+1;
    mi = m(ii);
    mj = m(jj);
    mij = m(ii)+m(jj);
    d00 = b300(k) + 1;
    d01 = b301(k) + 1;
    d10 = b310(k) + 1;
    d11 = b311(k) + 1;
    if length(beta) ==1
        dN = alpha*beta*N(ii)*N(jj);
    else
        dN = alpha*beta(k)*N(ii)*N(jj);
    end
    if dN > 0
        dMaggr(ii) = dMaggr(ii)-dN*mi;
        dMaggr(jj) = dMaggr(jj)-dN*mj;
        dMaggr(d00) = dMaggr(d00) + f00(k)*dN*mij; 
        dMaggr(d01) = dMaggr(d01) + f01(k)*dN*mij; 
        dMaggr(d10) = dMaggr(d10) + f10(k)*dN*mij; 
        dMaggr(d11) = dMaggr(d11) + f11(k)*dN*mij;
    end
end


%% Degradation

dNrem = remin(:).*N(:);
dNrom = -dNrem + [dNrem(2:end); 0]; % the last part corresponds to a shift of 1 to have the Z+1 * N(Z+1) 
dNrom(Nd:Nd:end) = -dNrem(Nd:Nd:end); % fix for bins: 10,20,30... (high density)--> nothing comes in
dMremin = dNrom(:).*m(:); % go back to dM (mass multiplying by mass per particle)
dMremin(1:Nd:end) = dMremin(1:Nd:end)- Rrate*M(1:Nd:end);  % ensure remin make mass to leave the system at Z=1 (low density)


%% Fragmentation
dMfrag = - pfrag(:).*M(:);
% dMfrag2 = - pfrag(:).*reshape(M,Nd,Nr);

dMfrag(1:Nd) = 0;

for i  = 1:Nr-1
    io = 1 + (i-1)*Nd:i*Nd; 
    for j = i:Nr-1
        jo = 1 + j*Nd:(j+1)*Nd;
        dMfrag(io) = dMfrag(io) - dMfrag(jo)/j;
    end
end


%% Sinking
dMsink = - M.*w(:)/H;


%% Collecting
dMdt = prod(:) + dMaggr + dMfrag + dMremin + dMsink;

end





