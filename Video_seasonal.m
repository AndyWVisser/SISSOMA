v = VideoWriter('Fig1_tempDepend_Rrate20', 'MPEG-4');
v.FrameRate = 12;
open(v);

Mo = sim.Mo;
Nd = sim.Nd;
Nr = sim.Nr;
w = sim.w;
mdry = sim.Mdry;
x = sim.x;
z = sim.z;

Rrate = sim.Rrate
remin_grid = sim.remin_grid;

% t = 1095:length(Mo(:,1));

T_mat = [10 10 linspace(10,20,182) linspace(20,10,182) linspace(10,20,182) linspace(20,10,182) linspace(10,20,182) linspace(20,10,182) linspace(10,20,182) linspace(20,10,182) linspace(10,20,182) linspace(20,10,182) 10 10 10 ];

for i = 1095:length(Mo(:,1))
    Mdry = reshape(Mo(i,:),Nd,Nr);
    Flux = Mdry.*w; % mass/ area/ time
    BFlux = sum(Flux,1);
    N = Mdry./(mdry);

    remin_c = Rrate(i) .* remin_grid;
    % remin_c = Rrate .* remin_grid;
    prod = sim.prod*(1-cos(2*pi*i/365))/2;


    figure(1);
    subplot(3,2,1); imagesc(x,z,prod); colorbar;  title('productivity [µgC/m2/day]'); axis xy  ; clim([0 2000])
    subplot(3,2,2); imagesc(x,z,Mdry); colorbar;  title('dry mass [µgC/m3]'); axis xy  ; clim([0 8*10^4]);
    subplot(3,2,3); imagesc(x,z,Flux); colorbar; title('flux [µgC/m2/day]'); axis xy  ; clim([0 8000])
    subplot(3,2,4); imagesc(x,z,remin_c); colorbar; title('Rrate [day-1]'); axis xy  ; clim([0 0.55])
    subplot(3,2,5:6);

    yyaxis left
    plot(x+.5,BFlux,'o-'); ylabel('size integrated flux') ; ylim([0 18000])

    yyaxis right
    plot(x(end),T_mat(i),'o-'); ylabel('Temperature') ; ylim([10 20])

   

    % subplot(2,2,4); semilogy(t(2:end),RMSE,'.');  title('root mean error');
    %     frame = getframe(gcf);
    %     writeVideo(v, frame);
    sgtitle(i)
    frame = getframe(gcf);
    writeVideo(v, frame);
end
close(v)