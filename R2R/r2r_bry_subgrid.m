function limits = new_bry_subgrid(parentgrid,childgrid,obcflag);
%
%   Find lower and upper index in i and j for the minimal
%   parent grid that contains all of the boundary bnd of
%   the child grid.
%
%   (c) 2007,  Jeroen Molemaker
%
%--------------------------------------------------------------

% get topography data from parentgrid
lonp = ncread(parentgrid,'lon_rho')';
latp = ncread(parentgrid,'lat_rho')';
lonp = double(lonp)+360;
latp = double(latp);
[Mp,Lp]=size(lonp);
tri_par = delaunay(lonp,latp);
%tri_par = delaunayTriangulation(reshape(lonp,Mp*Lp,1),reshape(latp,Mp*Lp,1));

limits = zeros(4,4);
for bnd = 1:4
 if ~obcflag(bnd)
    continue
 end
 if bnd == 1 %% South
   lonc = ncread(childgrid,'lon_rho')';lonc = lonc(1:2,1:end);
   latc = ncread(childgrid,'lat_rho')';latc = latc(1:2,1:end);
 end
 if bnd == 2 %% East
   lonc = ncread(childgrid,'lon_rho')';lonc = lonc(1:end,end-1:end);
   latc = ncread(childgrid,'lat_rho')';latc = latc(1:end,end-1:end);
 end
 if bnd == 3 %% North
   lonc = ncread(childgrid,'lon_rho')';lonc = lonc(end-1:end,1:end);
   latc = ncread(childgrid,'lat_rho')';latc = latc(end-1:end,1:end);
 end
 if bnd == 4 %% West
   lonc = ncread(childgrid,'lon_rho')';lonc = lonc(1:end,1:2);
   latc = ncread(childgrid,'lat_rho')';latc = latc(1:end,1:2);
 end

 Xp = [lonp(:) latp(:)];
 Xc = [lonc(:) latc(:)];
 [Mc,Lc] = size(lonc);
 size([reshape(lonc,Mc*Lc,1) reshape(latc,Mc*Lc,1) ]) 
 size(Xc)
 size(tri_par)
 size(Xp)
 t   = tsearch(lonp(:),latp(:),tri_par,lonc(:),latc(:));
%[nyc,nxc]= size(lonc);
%t   = squeeze(pointLocation(tri_par,reshape(lonc,nxc*nyc,1),reshape(latc,nxc*nyc,1)));

% Fix to deal with child points that are outside parent grid (those points should be masked!)
% if 0
  if (length(t(~isfinite(t)))>0);
    disp('Warning in new_bry_subgrid: outside point(s) detected.');
    [lonc,latc] = fix_outside_child(lonc,latc,t);
    Xc = [lonc(:) latc(:)];
    [t,pn] = tsearchn(Xp,tri_par,Xc);
%   [Mp,Lp] = size(lonp);
%   [Mc,Lc] = size(lonc);
%   Xp    = [reshape(lonp,Mp*Lp,1) reshape(latp,Mp*Lp,1) ]; Xp = double(Xp);
%   Xc    = [reshape(lonc,Mc*Lc,1) reshape(latc,Mc*Lc,1) ]; Xc = double(Xc);
%   plot(Xp(:,1),Xp(:,2),'.k')
%   hold on
%   plot(Xc(:,1),Xc(:,2),'.r')
%   hold off
%   drawnow
%   size(Xp)
%   size(Xc)
%   error
%   [t,pn] = tsearchn(Xp,tri_par,Xc);
  end;

  index       = tri_par(t,:);
% index       = tri_par.Triangulation(t,:);
  [idxj,idxi] = ind2sub([Mp Lp], index);

  limits(bnd,1)  = min(min(idxi));
  limits(bnd,2)  = max(max(idxi));
  limits(bnd,3)  = min(min(idxj));
  limits(bnd,4)  = max(max(idxj));

  if 0   %% debugging
   i0 = limits(bnd,1);
   i1 = limits(bnd,2);
   j0 = limits(bnd,3);
   j1 = limits(bnd,4);
   lonps= lonp(j0:j1,i0:i1);
   latps= latp(j0:j1,i0:i1);

   figure(bnd)
   plot(lonps,latps,'.k')
   hold on
   plot(lonc,latc,'.r')
   hold off
   error('debugging in bry_subgrid')
 end
end

