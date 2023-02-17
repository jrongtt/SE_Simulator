% Simulates Fractional Probabiltiy spent in finite well veruss potential 
% of well
% Takes ~= 30 minutes to run

tmax = 0.10;
level = 9;
lambda = 0.01;
idtype = 1; 
idpar = [0.40, 0.075, 0.0];
vtype = 1; % (rectangular well)

x1 = 0.6;
x2 = 0.8;




N = 251;

lnV = linspace(2, 10, N);
V = -exp(lnV);
Fe_vect = zeros(1, length(lnV));

for i = 1:length(lnV)
    vpar = [0.6, 0.8, V(i)]
    [x, t, psi, psire, psiim, psimod, prob, v] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);

    [mi, xmin_index] = min(abs(x-0.6));
    x1 = x(xmin_index);

    [ma, xmax_index] = min(abs(x-0.8));
    x2 = x(xmax_index);


    P_x1 =  mean(prob(:,xmin_index)) / mean(prob(:,nx));
    P_x2 =  mean(prob(:,xmax_index)) / mean(prob(:,nx));

    Fe = (P_x2 - P_x1) / (x2 - x1)
    Fe_vect(i) = Fe;

end



f = figure;
hold on;


ln_Fe = log(Fe_vect);

plot(lnV, ln_Fe);


xlabel("ln(V)", 'FontSize', 16);
ylabel("ln(Fe) Fraction Probability Potential",  'FontSize', 16);

title("Well Survey for Different Potentials", 'FontSize', 18)


hold off;
