ro = 1;
rMax = 1E6;
Nr = 30;
delta = exp(log(rMax/ro)/(Nr-1));
r_mean = [0.795519545878019	1.28099491690990	2.06273747219859	3.32154782430842	5.34856233421151	8.61258682882350	13.8685215295260	22.3319535974070	35.9602968791574	57.9055005643254	93.2430287456417	150.145708523880	241.773933037237	389.319383623924	626.906220043034	1009.48328097706	1625.53259481499	2617.53341199318	4214.91465920465	6787.11702512720	10929.0368221747	17598.6129925632	28338.3782396651	45632.2144019919	73479.8220850470	118321.767295485	190529.048909821	306801.692606786	494031.115596132	795519.545878020];

x = -1:28;
rLOW = ro*delta.^x

xx = 0:29;
rUP = ro*delta.^xx

dr = rUP*(1-1/delta);

rLOW1 = rUP/delta


% for i = 1:length(r_mean)
% 
%     rLOW1(i) = rUP/delta
%     % rLOw(i) = r_mean(i)/delta;
% end

xx = ones(30,1);
figure
semilogx(r_mean,xx, 'x', 'LineWidth', 1.2)
hold on
semilogx(rUP,xx, 'o', 'LineWidth', 3.9)
semilogx(rLOW,xx, 'o', 'LineWidth', 1.5)
xlabel('radius [\mum]')
legend('r_{mean}', 'Upper limit', 'Lower limit')

% dr = r*(1-1/sim.p.delta);

%%
x=0:29;
rUP = ro*delta.^x;
rMEAN = r_mean;

xx = ones(30,1);
figure
plot(rUP,xx, 'x')
hold on
plot(rMEAN,xx,'o')

%% Same as NUM

deltax = (log(1E6) - log(1)) / (30-1)

for i = 1:30

    x = log(1) + (i-1)*deltax
    m(i) = exp(x);
    mLow(i) = exp(x - 0.5*deltax);
    mUp(i) = exp(x + 0.5*deltax);
    mDelta(i) = mUp(i)-mLow(i);

end

xx = ones(30,1);
figure
semilogx(m,xx, 'x')
hold on
semilogx(mUp,xx, 'o')
semilogx(mLow,xx, '_')



%%
x =-1:28

rMINE = 1*delta.^x*(delta-1)./log(delta)
