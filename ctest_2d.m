% Script that produces the 2-D convergence plot discussed
% in the write up 
% This script takes ~= 30-60 min to run


idtype = 1;
vtype = 0;
tmax = 0.05;
lambda = 0.05;
lmin = 6;
lmax = 9;
vpar = 0;
level = 6
idpar = [0, 0.5, 0.1, 0.1, 0.8, 0];

[x6 y6 t6 psi6 psire6 psiim6 psimod6 v6] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);

level = 7

[x7 y7 t7 psi7 psire7 psiim7 psimod7 v7] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);

psi7_halved = (psi7(1:2:end, 1:2:end, 1:2:end));

psi7_6 = (psi7_halved - psi6);



level = 8

[x8 y8 t8 psi8 psire8 psiim8 psimod8 v8] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);

psi8_halved = (psi8(1:2:end, 1:2:end, 1:2:end));

psi8_7 = (psi8_halved - psi7);

level = 9

[x9 y9 t9 psi9 psire9 psiim9 psimod9 v9] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);

psi9_halved = (psi9(1:2:end, 1:2:end, 1:2:end));

psi9_8 = (psi9_halved - psi8);


rms_values7_6 = zeros(1, 9);

for i =1:length(t6)
    s = (psi7_6(i,:,:));
    rms_values7_6(i) = rms(rms(s, 2));
end

rms_values8_7 = zeros(1, 9);

for i =1:length(t7)
    s = squeeze(psi8_7(i,:,:));
    rms_values8_7(i) = rms(rms(s, 2));
end


rms_values9_8 = zeros(1, 9);
for i =1:length(t8)
    s = squeeze(psi9_8(i,:,:));
    rms_values9_8(i) = rms(rms(s, 2));
end





h = figure;
hold on;

plot(t6, rms_values7_6,"red");
plot(t7, rms_values8_7,"blue");
plot(t8, rms_values9_8,"green");

hold off;
xlabel("Time", 'FontSize', 16);
ylabel("RMS Value of dpsi",  'FontSize', 16);
%text(x, y, 'Figure', 'FontSize', 12);
title("2D Convergence plot for l6,l7,l8 dpsi", 'FontSize', 18)


legend({'dpsi6','4* dpsi7', '16*dpsi8'},'Location','southwest', 'FontSize', 16)

hold off

