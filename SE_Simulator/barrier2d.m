%Script that animates a Boosted Gaussian Wavefunction that encounters a
% rectangular barrier in 2D

level = 7;
tmax = 0.02;
lambda = 0.02;
idtype = 1;
vtype = 1;
vpar = [0, 1, 0.7, 0.8, 400]; % Barrier DONT CHANGE
idpar = [0.5, 0.4, 0.15, 0.15, 0, -15]; % Barrier



[x y t psi psire psiim psimod v] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);

nt = length(t);
obj = VideoWriter('barrier2dAnimation.avi');
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
    caxis([-1,1])
    title('Rectangular Barrier','FontSize', 18)
    hbar = colorbar;
    ylabel(hbar,'Abs(Psi)', 'FontSize', 16);
    hold off;
    f = getframe(gcf);
    writeVideo(obj, f);
    
end

