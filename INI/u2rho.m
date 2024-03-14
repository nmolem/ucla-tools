function [ur] = u2rho(u)
%
%
%  (c) 2020, Jeroen Molemaker, UCLA
%

 if ndims(u)==2
   [nx,ny] =size(u);
   ur = zeros(nx+1,ny);
   ur(2:nx,:) = 0.5*(u(1:nx-1,:) + u(2:nx,:));
   ur(1,:)    = u(1,:);
   ur(nx+1,:) = u(nx,:);
 elseif ndims(u)==3
   [nx,ny,nz]=size(u);
   ur = zeros(nx+1,ny,nz);
   ur(2:nx,:,:) = 0.5*(u(1:nx-1,:,:) + u(2:nx,:,:));
   ur(1,:,:)    = u(1,:,:);
   ur(nx+1,:,:) = u(nx,:,:);
 else
   error('only for 2 or 3 dimensions')
 end

return

