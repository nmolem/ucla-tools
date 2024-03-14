function [vr] = v2rho(v)
%
%
%  2007, Jeroen Molemaker, UCLA
%
 if ndims(v)==2
  [nx,ny]=size(v);
  vr = zeros(nx,ny+1);

  vr(:,2:ny) =0.5*(v(:,1:ny-1)+v(:,2:ny));
  vr(:,  1) =v(:,1);
  vr(:,ny+1)=v(:,ny);
 else
  [nx,ny,nz]=size(v);
  vr = zeros(nx,ny+1,nz);

  vr(:,2:ny,:) =0.5*(v(:,1:ny-1,:)+v(:,2:ny,:));
  vr(:,  1,:) =v(:,1,:);
  vr(:,ny+1,:)=v(:,ny,:);

 end

return

