% Script that produces the three 1-D convergence plots discussed
% in the write up 


idtype = 0;
vtype = 0;
idpar = [3];
tmax = 0.25;
lambda = 0.1;
lmin = 6;
lmax = 0;
vpar = 0;

level = 6;

[x6, t6, psi6, psire6, psiim6, psimod6, prob6, v6] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);

level = 7;

[x7, t7, psi7, psire7, psiim7, psimod7, prob7, v7] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);

psi7_halved = (psi7(1:2:end, 1:2:end));


level = 8;

[x8, t8, psi8, psire8, psiim8, psimod8, prob8, v8] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);

psi8_halved = (psi8(1:2:end, 1:2:end));

level = 9;

[x9, t9, psi9, psire9, psiim9, psimod9, prob9, v9] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);

psi9_halved = (psi9(1:2:end, 1:2:end));


rms7_6 = zeros(length(t6), 1);
for i=1:length(t6)
    s = squeeze(psi7_halved(i,:) - psi6(i,:));
    rms7_6(i) = rms(s);

end

rms8_7 = zeros(length(t7), 1);
for i=1:length(t7)
    s = squeeze(psi8_halved(i,:) - psi7(i,:));
    rms8_7(i) = rms(s);

end

rms9_8 = zeros(length(t8), 1);
for i=1:length(t8)
    s = squeeze(psi9_halved(i,:) - psi8(i,:));
    rms9_8(i) = rms(s);

end




%plot(x6, v7_6);





h = figure;
hold on;

plot(t6, rms7_6,"red");
plot(t7, 4*rms8_7,"blue");
plot(t8, 16*rms9_8,"green");

hold off;
xlabel("Time", 'FontSize', 16);
ylabel("RMS Value of dpsi",  'FontSize', 16);
%text(x, y, 'Figure', 'FontSize', 12);
title("Convergence plot for l6,l7,l8 dpsi", 'FontSize', 18)


legend({'dpsi6','4* dpsi7', '16*dpsi8'},'Location','southwest', 'FontSize', 16)

hold off


t = figure;
hold on;




psi6_exact = sin(idpar(1)*pi*x6);
psi7_exact = sin(idpar(1)*pi*x7);
psi8_exact = sin(idpar(1)*pi*x8);
psi9_exact = sin(idpar(1)*pi*x9);



rms6 = zeros(length(t6), 1);
for i=1:length(t6)
    s = squeeze(psi6_exact - psi6(i,:));
    rms6(i) = rms(s);

end

rms7 = zeros(length(t7), 1);
for i=1:length(t7)
    s = squeeze(psi7_exact - psi7(i,:));
    rms7(i) = rms(s);

end

rms8 = zeros(length(t8), 1);
for i=1:length(t8)
    s = squeeze(psi8_exact - psi8(i,:));
    rms8(i) = rms(s);

end

rms9 = zeros(length(t9), 1);
for i=1:length(t9)
    s = squeeze(psi9_exact - psi9(i,:));
    rms9(i) = rms(s);

end

plot(t6, rms6,"red");
plot(t7, rms7,"blue");
plot(t8,  rms8,"green");
plot(t9,  rms9,"black");





xlabel("Time", 'FontSize', 16);
ylabel("RMS Value of ||E(psi)||",  'FontSize', 16);
%text(x, y, 'Figure', 'FontSize', 12);
title("Exact Convergence for ||E(psi)|| for levels: l6,l7,l8,l9", 'FontSize', 18)


legend({'Level 6','Level 7', 'Level 8', 'Level 9'},'Location','southwest', 'FontSize', 16)


idtype = 1;
vtype = 0;
idpar = [0.50 0.075 0.0];
tmax = 0.01;
lambda = 0.01;
lmin= 6;
lmax= 9;
vpar = [];


level = 6;

[x6, t6, psi6, psire6, psiim6, psimod6, prob6, v6] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);

level = 7;

[x7, t7, psi7, psire7, psiim7, psimod7, prob7, v7] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);

psi7_halved = (psi7(1:2:end, 1:2:end));


level = 8;

[x8, t8, psi8, psire8, psiim8, psimod8, prob8, v8] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);

psi8_halved = (psi8(1:2:end, 1:2:end));

level = 9;

[x9, t9, psi9, psire9, psiim9, psimod9, prob9, v9] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);

psi9_halved = (psi9(1:2:end, 1:2:end));




psi7_6 = (psi7_halved - psi6);
psi8_7 = (psi8_halved - psi7);
psi9_8 = (psi9_halved - psi8);



rms7_6 = zeros(length(t6), 1);
for i=1:length(t6)
    s = squeeze(psi7_halved(i,:) - psi6(i,:));
    rms7_6(i) = rms(s);

end

rms8_7 = zeros(length(t7), 1);
for i=1:length(t7)
    s = squeeze(psi8_halved(i,:) - psi7(i,:));
    rms8_7(i) = rms(s);

end

rms9_8 = zeros(length(t8), 1);
for i=1:length(t8)
    s = squeeze(psi9_halved(i,:) - psi8(i,:));
    rms9_8(i) = rms(s);

end




ss = figure;
hold on;

plot(t6, rms7_6,"red");
plot(t7, 4*rms8_7,"blue");
plot(t8, 16*rms9_8,"green");

xlabel("Time", 'FontSize', 16);
ylabel("RMS Value of dpsi",  'FontSize', 16);
title("Convergence plot for l6,l7,l8 dpsi", 'FontSize', 18)


legend({'dpsi6','4* dpsi7', '16*dpsi8'},'Location','southwest', 'FontSize', 16)




hold off;
xlabel("Time")
ylabel("RMS for dpsi")
title("Second convergence for l6,l7,l8 dpsi")

