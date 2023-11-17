tt = t(end-730:end);

for i = 1:length(tt)

    prod = Ptotal/Nd/H .* (1-cos(2*pi*tt(i)/365))/2;
    Mdry2 = squeeze(Mo(tt(i),:));
    [dMdto,dMsink,dMremin,dMfrag,dMaggreg] = interaxseason(tt(i),Mdry2(:),mdry,bi,bj,Nr,Nd,b300,b301,b310,b311,f00,f01,f10,f11,alpha,beta,w,H,prod,Rrate,pfrag, remin_grid);
    
    adMdt(i,:) = dMdto;
    adMsink(i,:) = dMsink;
    adMremin(i,:) = dMremin;
    aprod(i,:) = Ptotal/Nd/H .* (1-cos(2*pi*tt(i)/365))/2;
    
end



%%


adMremin2 = sum(adMremin,2);
adMsink2 = sum(adMsink,2);

ares = aprod - (adMremin2 + adMsink2);

