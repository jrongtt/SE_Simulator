%Script that computs the fractional probability for particles exsisting to
%the right of barriers of varying potential
% This script takes ~= 30 minutes to run

tmax = 0.10;
level = 9;
nx = 2^level + 1;
lambda = 0.01;
idtype = 1; % Boosted Gaussian
idpar =  [0.40, 0.075, 20.0];
vtype = 1; %Rectangular Barrier

x1 = 0.8;
x2 = 1.0;

N = 251;

lnV = linspace(-2, 5, N);
V = exp(lnV);
Fe_vect = zeros(1, length(lnV));

for i = 1:length(lnV)
    i
    vpar = [0.6, 0.8, V(i)]
    [x, t, psi, psire, psiim, psimod, prob, v] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);
    
    [mi, xmin_index] = min(abs(x-x1));
    x1 = x(xmin_index);

    [ma, xmax_index] = min(abs(x-x2));
    x2 = x(xmax_index);
    
    P_x1 =  mean(prob(:,xmin_index)) / mean(prob(:,nx));
    P_x2 =  mean(prob(:,xmax_index)) / mean(prob(:,nx));

    Fe = (P_x2 - P_x1) / (x2 - x1);
    Fe_vect(i) = Fe;

end


f = figure;
hold on;
ln_Fe = log(Fe_vect);

plot(lnV, ln_Fe);


xlabel("ln(V)", 'FontSize', 16);
ylabel("ln(Fe) Fraction Probability Potential",  'FontSize', 16);
title("Barrier Survey for Different Potentials", 'FontSize', 18)


hold off;