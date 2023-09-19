
G = []; Nx = []; D=[]; Wx = [];
paramset = 'a'
switch paramset
    case {'a','selfsimilarity' }
        load sima1y0.mat; % self similarity
        legendtext = {'sima17','sima18','sima19','sima20','sima21','sima22','sima23'};
        legendtext0 = {'1.7','1.8','1.9','2.0','2.1','2.2','2.3'};
    case {'p','production'}
        load simp1year0.mat; % self similarity
        legendtext = {'simp01','simp02','simp05'};
        legendtext0 = {'1E5','2E5','5E5'};
    case {'s','stickiness'}
        load sims1y0.mat; % self similarity
        legendtext = {'sims01','sims02','sims03','sims04','sims05','sims06','sims07','sims08'}
        legendtext0 = {'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'}
    case {'e','dissipation'}
        load sime1y0.mat; % self similarity
        legendtext = {'sime02','sime03','sime04','sime05','sime06','sime07','sime08'}
        legendtext0 = {'1E-2','1E-3','1E-4','1E-5','1E-6','1E-7','1E-8'};
    otherwise
        return
end

sim = evalin('base',cell2mat(legendtext(1)));
x = sim.x;
z = sim.z;
r = sim.p.rMin*sim.p.delta.^x;
dr = r*(1-1/sim.p.delta);
Flux = @(sim) sum((sim.w).*(sim.M))/sim.p.H/sum(sum(sim.prod));
Nspec = @(sim) sum(sim.M./sim.m);
Dmean = @(sim) 1E9*sum(sim.N.*sim.d)./sum(sim.N);
Wmean = @(sim) sum(sim.N.*sim.w)./sum(sim.N);
G(1,:) = Flux(sim); Nx(1,:) = Nspec(sim)./dr; D(1,:) = Dmean(sim); Wx(1,:) = Wmean(sim);
for i = 2:length(legendtext)
    sim = evalin('base',cell2mat(legendtext(i)));
    G(i,:) = Flux(sim); Nx(i,:) = Nspec(sim)./dr; D(i,:) = Dmean(sim); Wx(i,:) = Wmean(sim);
end
sz = size(G);
out.g = G;
out.n = Nx;
out.d = D;
out.r = r;
out.l = legendtext; out.l0 = legendtext0;
switch paramset
    case {'a','selfsimilarity' }
        outa = out;
    case {'p','production'}
        outp = out;
    case {'s','stickiness'}
        outs = out;
    case {'e','dissipation'}
        oute = out;
    otherwise
        return
end

figure(1)
clf;
hg = semilogx(r,G,'LineWidth',1);
ag = gca; ag.FontSize = 12;
legend (legendtext0)

figure(2)
clf;
s = 1E13*r.^(-4); 
hn = loglog(r,Nx,r,s,'--','LineWidth',1);
an = gca; an.FontSize = 12;
set(hn(end),'Color',[0.7,0.7,0.7])
legend (legendtext0)
ylim([1E-10 1E15])
xlim([1E0 1E6])

figure(3)
clf;
s = 1E13*r.^(-4); 
hd = loglog(r,D,'LineWidth',1);
ad = gca; ad.FontSize = 12;
% set(hd(end),'Color',[0.7,0.7,0.7])
legend (legendtext0)
ylim([1E-2 1E2])
xlim([1E0 1E6])

figure(4)
clf;
rMesh = ones(sz(1),1)*r;

rfmax = sum(G.*rMesh,2)./sum(G,2);
wfmax = sum(G.*Wx,2)./sum(G,2);

hd = loglog(r,Wx,'-',rfmax,wfmax,'o');
ad = gca; ad.FontSize = 12;
% set(hd(end),'Color',[0.7,0.7,0.7])
legend (legendtext0)
ylim([1E-2 1E3])
xlim([1E0 1E6]);

figure(5)
clf;
s = 1E13*r.^(-4); 
hd = plot(1:sz(1),sum(G.*Wx,2)./sum(G,2),'o');
ad = gca; ad.FontSize = 12;
% set(hd(end),'Color',[0.7,0.7,0.7])
legend (legendtext0)
ylim([0 100])
%xlim([1E0 1E6]);