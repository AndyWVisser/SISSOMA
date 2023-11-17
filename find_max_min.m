maxM = -inf;
minM = inf;


for i = 1095:1825


    remin_c = remin(i) .* remin_grid;
    % remin_c = sim.prod*(1-cos(2*pi*i/365))/2;

    if max(max(remin_c)) > maxM
        maxM = max(max(remin_c));
    end

    if min(min(remin_c)) < minM
        minM = min(min(remin_c));
    end

end