function [dMdt,dMsink,dMremin,dMfrag] = interax(t,M,m,xzb,bi,bj,Nr,Nd,q,a,b300,b301,b310,b311,f00,f01,f10,f11,alpha,beta,w,L,H,prod,remin,rfact,pfrag,phi)
% returns the time rate of change of dry mass per size, excess-density interval. 
% 
%  
M(M<0) = 0;
m=m(:); 
N = M(:)./m(:);
dM = zeros(size(M));
dMremin = zeros(size(dM));
dMfrag = zeros(size(dM));

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
        dM(ii) = dM(ii)-dN*mi;
        dM(jj) = dM(jj)-dN*mj;
        dM(d00) = dM(d00) + f00(k)*dN*mij; 
        dM(d01) = dM(d01) + f01(k)*dN*mij; 
        dM(d10) = dM(d10) + f10(k)*dN*mij; 
        dM(d11) = dM(d11) + f11(k)*dN*mij;
    end
end
%% Degradation
dNrem = remin(:).*N(:);
dNrom = -dNrem + [dNrem(2:end); 0];
dNrom(Nd:Nd:end) = -dNrem(Nd:Nd:end);
dMremin = dNrom.*m(:);
% dMremin =  - rfact*M(:);
%% Fragmentation
% dMfrag = - pfrag(:).*M(:)./m(:); dMfrag(1:Nd) = 0;
% for i  = 1:Nr-1
%     io = 1 + (i-1)*Nd:i*Nd; 
%     for j = i:Nr-1
%         jo = 1 + j*Nd:(j+1)*Nd;
%         dMfrag(io) = dMfrag(io) - dMfrag(jo)/j;
%     end
% end
%% Fragmentation
dMfrag = - pfrag(:).*M(:); dMfrag(1:Nd) = 0;
for i  = 1:Nr-1
    io = 1 + (i-1)*Nd:i*Nd; 
    for j = i:Nr-1
        jo = 1 + j*Nd:(j+1)*Nd;
        dMfrag(io) = dMfrag(io) - dMfrag(jo)/j;
    end
end
%% Collecting
dMsink = - M.*w(:)/H;
dMdt = dM + prod(:) + dMsink + dMremin + dMfrag;
end