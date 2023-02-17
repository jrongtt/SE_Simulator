    function [x, t, psi, psire, psiim, psimod, prob, v] = ...
sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar)
% Inputs
%
% tmax: Maximum integration time
% level: Discretization level
% lambda: dt/dx
% idtype: Selects initial data type
% idpar: Vector of initial data parameters
% vtype: Selects potential type
% vpar: Vector of potential parametersm
%
% Outputs
%
% x: Vector of x coordinates [nx]
% t: Vector of t coordinates [nt]
% psi: Array of computed psi values [nt x nx]
% psire Array of computed psi_re values [nt x nx]
% psiim Array of computed psi_im values [nt x nx]
% psimod Array of computed sqrt(psi psi*) values [nt x nx]
% prob Array of computed running integral values [nt x nx]
% v Array of potential values [nx]


 % Define mesh and derived parameters ...

   nx = 2^level + 1;
   x = linspace(0.0, 1.0, nx);
   dx = x(2) - x(1);
   dt = lambda * dx;
   nt = round(tmax / dt) + 1;
   t = [0 : nt-1] * dt;
   v = zeros(1,nx);

   
    % Initialize solution, and set initial data ...
   psi = zeros(nt, nx);
   if idtype == 0
      psi(1, :) = sin(idpar(1)*pi*x);
   elseif idtype == 1
      
      psi(1, :) = exp(1i*idpar(3)*x) .* exp(-((x-idpar(1))/idpar(2)).^2);
   else
      fprintf('diff_1d_cn: Invalid idtype=%d\n', idtype);
      return
   end

   if vtype == 1

       [mi, xmin_index] = min(abs(x-vpar(1)));
       x1 = x(xmin_index);

       [ma, xmax_index] = min(abs(x-vpar(2)));
        x2 = x(xmax_index);

       for i = xmin_index : xmax_index
           v(i) = vpar(3);
       end
   end


   % Initialize storage for sparse matrix and RHS ...
   dl = zeros(nx,1);
   d  = zeros(nx,1);
   du = zeros(nx,1);
   f  = zeros(nx,1);

   % Set up tridiagonal system ...
   dl = 0.5 / dx^2 * ones(nx, 1);
   d = (((1i / dt) - (1 / dx^2)) * ones(nx,1)) - 0.5*v.';


  
  % d  = (1.0 / dt + 1.0 / dx^2) * ones(nx,1);

   du = dl;
   % Fix up boundary cases ...
   d(1) = 1.0;
   du(2) = 0.0;
   dl(nx-1) = 0.0;
   d(nx) = 1.0;
   % Define sparse matrix ...
   A = spdiags([dl d du], -1:1, nx, nx);
   prob = zeros(nt,nx);
   psimod = psi;
 

   for n = 1 : nt-1
      % Define RHS of linear system ...
   

      f(2:nx-1) = -0.5*( psi(n, 1:nx-2) - 2 * psi(n, 2:nx-1) + psi(n, 3:nx))./dx^2 + (1i* psi(n, 2:nx-1)./dt) + 0.5*v(2:nx-1).*psi(n, 2:nx-1);
      

      f(1) = 0.0;
      f(nx) = 0.0;
      % Solve system, thus updating approximation to next time 
      % step ...
      psi(n+1, :) = A \ f;
      psimod(n+1,:) = sqrt(psi(n+1, :).* conj(psi(n+1, :)));

      % Trapezoidal Integral   
      sum = 0;
      for i =1:nx-1
          sum = sum + (psi(n+1, i).*conj(psi(n+1, i)) + psi(n+1, i+1).*conj(psi(n+1, i+1))) * (x(i+1) - x(i));
          %sum = sum + ((abs(psi(n+1, i))+ abs(psi(n+1, i+1))) * (x(i+1) - x(i)));
          prob(n+1, i+1) = 0.5*sum;
    
      end
      

   end   
  
   

  
   
      

   %From solution u, get necessary outputs:
   
  
   psire = real(psi);
   psiim = imag(psi);




end