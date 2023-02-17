% Inputs
%
% tmax: Maximum integration time
% level: Discretization level
% lambda: dt/dx
% idtype: Selects initial data type
% idpar: Vector of initial data parameters
% vtype: Selects potential type
% vpar: Vector of potential parameters
%
% Outputs
%
% x: Vector of x coordinates [nx]
% y: Vector of y coordinates [ny]
% t: Vector of t coordinates [nt]
% psi: Array of computed psi values [nt x nx x ny]
% psire Array of computed psi_re values [nt x nx x ny]
% psiim Array of computed psi_im values [nt x nx x ny]
% psimod Array of computed sqrt(psi psi*) values [nt x nx x ny]
% v Array of potential values [nx x ny]

function [x y t psi psire psiim psimod v] = ...
    sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar)


    % Define mesh and derived parameters ...

    nx = 2^level + 1;
    ny = 2^level + 1;
    x = linspace(0.0, 1.0, nx);
    y = linspace(0.0, 1.0, ny);
    dx = x(2) - x(1);
    dy = y(2) - y(1);
    dt = lambda * dx;
    nt = round(tmax / dt) + 1;
    t = [0 : nt-1] * dt;
    v = zeros(nx, ny);

      % Initialize solution, and set initial data ...
   psi = zeros(nt,nx,ny);
   if idtype == 0
      for i = 1:nx
          for j = i:ny
              psi(1,i,j)= sin(idpar(1)*pi*x(i)) * sin(idpar(2)*pi*y(j));
              psi(1,j,i)= sin(idpar(1)*pi*x(j)) * sin(idpar(2)*pi*y(i));
          end
      end
   elseif idtype == 1
       for i = 1:nx
          for j = i:ny
              psi(1,j,i)= exp(-(1i*idpar(5)*x(i)))*exp(-(1i*idpar(6)*y(j))) * exp(-((x(i)-idpar(1))^2 / idpar(3)^2) + -((y(j)-idpar(2))^2 / idpar(4)^2));
              psi(1,i,j)= exp(-(1i*idpar(5)*x(j)))*exp(-(1i*idpar(6)*y(i))) * exp(-((x(j)-idpar(1))^2 / idpar(3)^2) + -((y(i)-idpar(2))^2 / idpar(4)^2));
          end
      end
       

   else
      fprintf('diff_1d_cn: Invalid idtype=%d\n', idtype);
      return
   end


    if vtype == 1
       xmin_index = find((x-vpar(1))>=0, 1, 'first');
       xmax_index = find((x-vpar(2))>=0, 1, 'first');
       ymin_index = find((y-vpar(3))>=0, 1, 'first');
       ymax_index = find((y-vpar(4))>=0, 1, 'first');
       

       for i = xmin_index : xmax_index
           for j = ymin_index : ymax_index
               v(j, i) = vpar(5);
           end
       end
    end

    if vtype == 2
        %j_p denotes the j prime in the project hand-out
        v = zeros(nx, ny);
        j_p = round((ny - 1)/4) + 1;
        for  i=1:nx
            if ((vpar(1) <= x(i) && x(i) <= vpar(2))) || ((vpar(3) <= x(i)) && (x(i) <= vpar(4)))
%                 v(i, j_p) = 0;
%                 v(i, j_p +1) = 0;

                v(j_p, i) = 0;
                v(j_p + 1, i) = 0;
            else
                v(j_p, i) = vpar(5);
                v(j_p + 1, i) = vpar(5);
                
            end
        end        
    end

    % Initialize storage for sparse matrix and RHS ...
 

   % Set up tridiagonal system ...
   dl = (-1i*dt)/ (2*dx^2) * ones(nx, 1);
   d = (1 + (1i*dt / dx^2)) * ones(nx,1);
   du = dl;
   % Fix up boundary cases ...
   d(1) = 1.0;
   du(2) = 0.0;
   dl(ny-1) = 0.0;
   d(ny) = 1.0;
   % Define sparse matrix ...
   A = spdiags([dl d du], -1:1, ny, ny);


% 
   dl2 = (-1i*dt)/ (2*dy^2) * ones(ny, 1);
   
   du2 = dl2;
   % Fix up boundary cases ...
   d2(1) = 1.0;
   du2(2) = 0.0;
   dl2(nx-1) = 0.0;
   d2(nx) = 1.0;
   % Define sparse matrix ...
   
   
   

   for n = 1: nt-1
       psi_half = zeros(nx, ny);       
       temp_f = zeros(nx, ny);

       f_right = squeeze(psi(n, 2:nx-1,2:ny-1)) + ...
           (1i*dt)/(2*dy^2)*(squeeze(psi(n,2:nx-1,3:ny))...
           -2*squeeze(psi(n,2:nx-1,2:ny-1)) + squeeze(psi(n,2:nx-1,1:ny-2))) ...
           - (1i*dt*0.5*v(2:nx-1,2:ny-1).*squeeze(psi(n,2:nx-1,2:ny-1)));
 

       temp_f(2:nx-1, 2:ny-1) = f_right;

       f_left = zeros(nx, ny);
       f_left(2:nx-1, 2:ny-1) = temp_f(2:nx-1, 2:ny-1) + (1i*dt*0.5/dx^2)*...
           (temp_f(3:nx, 2:ny-1) -2*temp_f(2:nx-1, 2:ny-1) + temp_f(1:nx-2, 2:ny-1));

       f = zeros(ny, 1);
       
       
       % Iterate through js
        for j = 2: ny-1
            f(2:nx-1) =  (f_left(2:nx-1, j));
            f(1) = 0;
            f(ny) = 0;

            
            psi_half(:,j) = A \ f;
            psi_half(nx, :) = 0;
            psi_half(1, :) = 0;
            psi_half(:, ny) = 0;
            psi_half(:, 1) = 0;
        end
        
          f2 = zeros(nx, 1);

          % Set up second tridiagonal system ...
          d2 = (1+ (1i*dt / dy^2)) * ones(ny,1)  + (1i*dt*0.5*v(:,j)) ;
          A2 = spdiags([dl2 d2 du2], -1:1, nx, nx);


        for i = 2:nx-1
            f2(2:ny-1) = psi_half(i, 2:ny-1);
            f2(1) = 0;
            f2(nx) = 0;
            psi(n+1, i, :) = A2 \ f2;

        end

   end

   %From solution psi, get necessary outputs:
   psire = real(psi);
   psiim = imag(psi);
   psimod = psi .* (psire -psiim);
  


end








