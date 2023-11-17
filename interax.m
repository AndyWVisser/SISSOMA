function [dMdt,dMsink,dMremin,dMfrag] = interax(t,M,mdry,bi,bj,Nr,Nd,b300,b301,b310,b311,f00,f01,f10,f11,alpha,beta,w,H,prod,remin,rfact,pfrag)
% returns the time rate of change of dry mass per size, excess-density interval. 

% Dry mass [µgC/m3]
M(M<0) = 0; 

% Dry mass per particle reference --> mdry

% Number of particles per volume [#/ m3]
N = M(:)./mdry(:);

% Pre-allocate memory for the aggregation matrix
dMaggr = zeros(size(M));


%% Aggregation
for k = 1:length(bi)
    ii = bi(k)+1;
    jj = bj(k)+1;
    mi = mdry(ii);
    mj = mdry(jj);
    mij = mdry(ii)+mdry(jj);
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


%% Remineralization
dNrem = remin(:).*N(:);
dNrom = -dNrem + [dNrem(2:end); 0];
dNrom(Nd:Nd:end) = -dNrem(Nd:Nd:end);
dMremin = dNrom.*mdry(:);


%% Fragmentation
dMfrag = - pfrag(:).*M(:); 
dMfrag(1:Nd) = 0;

for i  = 1:Nr-1
    io = 1 + (i-1)*Nd:i*Nd; 
    for j = i:Nr-1
        jo = 1 + j*Nd:(j+1)*Nd;
        dMfrag(io) = dMfrag(io) - dMfrag(jo)/j;
    end
end


%% Sinking
dMsink = -M.*w(:)/H;


%% Collecting
dMdt = prod(:) + dMaggr + dMfrag + dMremin + dMsink;



end