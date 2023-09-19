outa = collectplot(1);
outs = collectplot(2);
oute = collectplot(3);
outf = collectplot(4);

figure(1)
clf;
tfigs = tiledlayout(4,2)
nexttile
out = outa;
hg = semilogx(out.r,out.g,'LineWidth',2);
ag = gca; ag.FontSize = 12; ag.LineWidth = 1;
ag.YAxis.TickLength = [0.02 0.02]
ag.XAxis.TickLength = [0.02 0.02]

nexttile
s3 = 1E11*out.r.^(-3);
s4 = 1E10*out.r.^(-4);
hn = loglog(out.r,out.n,out.r,s3,'--',outa.r,s4,'--','LineWidth',2);
an = gca; an.FontSize = 12; an.LineWidth = 1
an.YAxisLocation = 'right'
an.YAxis.TickLength = [0.02 0.02]
an.XAxis.TickLength = [0.02 0.02]
set(hn(end),'Color',[0.7,0.7,0.7])
set(hn(end-1),'Color',[0.7,0.7,0.7])
legend (out.leg)
ylim([1E-10 1E10])
xlim([1E0 1E6])
tfigs.TileSpacing = 'compact';
tfigs.Padding = 'compact';

nexttile
out = outs;
hg = semilogx(out.r,out.g,'LineWidth',2);
ag = gca; ag.FontSize = 12; ag.LineWidth = 1;
ag.YAxis.TickLength = [0.02 0.02]
ag.XAxis.TickLength = [0.02 0.02]

nexttile
s3 = 1E11*out.r.^(-3);
s4 = 1E10*out.r.^(-4);
hn = loglog(out.r,out.n,out.r,s3,'--',outa.r,s4,'--','LineWidth',2);
an = gca; an.FontSize = 12; an.LineWidth = 1
an.YAxisLocation = 'right'
an.YAxis.TickLength = [0.02 0.02]
an.XAxis.TickLength = [0.02 0.02]
set(hn(end),'Color',[0.7,0.7,0.7])
set(hn(end-1),'Color',[0.7,0.7,0.7])
legend (out.leg)
ylim([1E-10 1E10])
xlim([1E0 1E6])
tfigs.TileSpacing = 'compact';
tfigs.Padding = 'compact';

nexttile
out = oute;
hg = semilogx(out.r,out.g,'LineWidth',2);
ag = gca; ag.FontSize = 12; ag.LineWidth = 1;
ag.YAxis.TickLength = [0.02 0.02]
ag.XAxis.TickLength = [0.02 0.02]

nexttile
s3 = 1E11*out.r.^(-3);
s4 = 1E10*out.r.^(-4);
hn = loglog(out.r,out.n,out.r,s3,'--',outa.r,s4,'--','LineWidth',2);
an = gca; an.FontSize = 12; an.LineWidth = 1
an.YAxisLocation = 'right'
an.YAxis.TickLength = [0.02 0.02]
an.XAxis.TickLength = [0.02 0.02]
set(hn(end),'Color',[0.7,0.7,0.7])
set(hn(end-1),'Color',[0.7,0.7,0.7])
legend (out.leg)
ylim([1E-10 1E10])
xlim([1E0 1E6])
tfigs.TileSpacing = 'compact';
tfigs.Padding = 'compact';

nexttile
out = outf;
hg = semilogx(out.r,out.g,'LineWidth',2);
ag = gca; ag.FontSize = 12; ag.LineWidth = 1;
ag.YAxis.TickLength = [0.02 0.02]
ag.XAxis.TickLength = [0.02 0.02]

nexttile
s3 = 1E11*out.r.^(-3);
s4 = 1E10*out.r.^(-4);
hn = loglog(out.r,out.n,out.r,s3,'--',outa.r,s4,'--','LineWidth',2);
an = gca; an.FontSize = 12; an.LineWidth = 1
an.YAxisLocation = 'right'
an.YAxis.TickLength = [0.02 0.02]
an.XAxis.TickLength = [0.02 0.02]
set(hn(end),'Color',[0.7,0.7,0.7])
set(hn(end-1),'Color',[0.7,0.7,0.7])
legend (out.leg)
ylim([1E-10 1E10])
xlim([1E0 1E6])
tfigs.TileSpacing = 'compact';
tfigs.Padding = 'compact';