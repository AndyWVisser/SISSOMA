figure(1)
clf;
sim = sima20;
x = sim.x;
r = sim.p.rMin*sim.p.delta.^x;
z = sim.z;
sz = size(sim.M);
t = sim.t;
prod = sum(sum(sim.prod))/2;

dM = sim.dM(end-365:end,:);
w = sim.w;
m = sim.m;
H = sim.p.H;
F = 2*w.*reshape(dM(1,:),sz)/H/prod;
hF = surf(log10(r),z,F);
time = 0;
zlim([0 .5]);
xlabel('aggregate size')
ylabel('excess density')
zlabel('flux')
grid off
hT = text(1,30,0.4,['time = ' num2str(time,'%.2f')]);

dM = sim.dM;


for i = 1:length(t)
    time = i/365;
    F = 2*w.*reshape(dM(i,:),sz)/H/prod;
    set(hF,'zdata',F);
    set(hT,'String',['time = ' num2str(time,'%.2f')]);
    drawnow
end


   


