x = 1:1825; % SOS!!
y = Rrate;
Rrate2 = interp1(x, y, 185.956)

% remin = Rrate2 .* remin_grid;

%%

Temp = 10;
tt = 1:730;
Ts = Temp*(3-cos(2*pi*tt/365))/3

figure
plot(Ts)

%%

T0 = 10;
tt = 1:Tmax;
T_input = T0*(3-cos(2*pi*tt/365))/3;
figure
plot(T_input)

%%

T1 = 6:0.1:13;
Rrate_ref10 = 0.07; % remineralisation rate (1/day) (Serra-Pompei (2022)) @10 degrees
Q10 = 2;
Rrate1 = Rrate_ref10 .* fTemp(Q10, T1);


%%

x = 1:71; % SOS!!
y = Rrate1;


n = 50;
R = [1 71];
z = rand(n,1)*range(R)+min(R)

for i = 1:length(z)

    Rrate2(i) = interp1(x, y, z(i));

end


figure
plot(Rrate1,'o')
hold on
plot(Rrate2,'x')
hold on



