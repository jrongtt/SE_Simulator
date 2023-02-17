%Script that creates a Double slit and intilaizes a particle with y
%momentum that travels through slit (this takes a minute or two to run)

level = 8;
tmax = 0.007;
lambda = 0.02;
idtype = 1;
vtype = 2;
vpar = [0.38, 0.41, 0.59, 0.62, 500]; 

idpar = [0.5, 0.05, 0.07, 0.07, -2, -10];
[x y t psi psire psiim psimod v] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);
nt = length(t);
obj = VideoWriter('doubleSlit2dAnimation.avi');
obj.Quality = 100;
obj.FrameRate = 30;
open(obj)

figure('Position', [200 200 600 600])

for n=1:nt
    contourf(x,y,abs(squeeze(psi(n,:,:))), 20);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    zlim([-1, 1]);
    caxis([-1,1]);
    title('Double Slit','FontSize', 18);
    hbar = colorbar;
    ylabel(hbar,'Abs(Psi)', 'FontSize', 16);
    hold off;
    f = getframe(gcf);
    writeVideo(obj, f);
end

