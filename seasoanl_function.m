t = 1:730


pp = (1-cos(2*pi*t/365))/2.5 + 0.2
min(pp)
plot(pp)
hold on