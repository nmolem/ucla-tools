function limits = new_bry_subgrid(parentgrid,childgrid,obcflag);
%
%   Find lower and upper index in i and j for the minimal
%   parent grid that contains all of the boundary bnd of
%   the child grid.
%
%   (c) 2007,2014  Jeroen Molemaker
%
%--------------------------------------------------------------

% get topography data from parentgrid
nc=netcdf(parentgrid,'nowrite');
lonp = nc{'lon_rho'}(:);
latp = nc{'lat_rho'}(:);
close(nc);
[Mp,Lp]=size(lonp);
!!!!!!!!!!!!!!!!!!!        check x vs. y definition
%tri_par = delaunay(lonp,latp);
tri_par = DelaunayTri(reshape(lonp,Mp*Lp,1),reshape(latp,Mp*Lp,1));

limits = zeros(4,4);
for bnd = 1:4
 if ~obcflag(bnd)
    continue
 end
 nc = netcdf(childgrid,'nowrite');
 if bnd == 1 %% South
   lonc = ncread(childgrid,'lon_rho',[1 1],[inf 2]);
   latc = ncread(childgrid,'lat_rho',[1 1],[inf 2]);
%  lonc = nc{'lon_rho'}(1:2,1:end);
%  latc = nc{'lat_rho'}(1:2,1:end);
 end
 if bnd == 2 %% East
   lonc = ncread(childgrid,'lon_rho',[L 1],[2 inf]);
   latc = ncread(childgrid,'lat_rho',[L 1],[2 inf]);
%  lonc = nc{'lon_rho'}(1:end,end-1:end); 
%  latc = nc{'lat_rho'}(1:end,end-1:end);
 end
 if bnd == 3 %% North
   lonc = nc{'lon_rho'}(end-1:end,1:end); 
   latc = nc{'lat_rho'}(end-1:end,1:end);
 end
 if bnd == 4 %% West
   lonc = nc{'lon_rho'}(1:end,1:2); 
   latc = nc{'lat_rho'}(1:end,1:2);
 end
 close(nc);

%t   = squeeze(tsearch(lonp,latp,tri_par,lonc,latc));
 [nyc,nxc]= size(lonc);
 t   = squeeze(pointLocation(tri_par,reshape(lonc,nxc*nyc,1),reshape(latc,nxc*nyc,1)));

% Fix to deal with child points that are outside parent grid (those points should be masked!)
  if (length(t(~isfinite(t)))>0);
    disp('Warning in new_bry_subgrid: outside point(s) detected.');
    [lonc,latc] = fix_outside_child(lonc,latc,t);
    t  = squeeze(tsearch(lonp,latp,tri_par,lonc,latc));
  end;

% index       = tri_par(t,:);
  index       = tri_par.Triangulation(t,:);
  [idxj,idxi] = ind2sub([Mp Lp], index);

  limits(bnd,1)  = min(min(idxi));
  limits(bnd,2)  = max(max(idxi));
  limits(bnd,3)  = min(min(idxj));
  limits(bnd,4)  = max(max(idxj));

  if 1   %% debugging
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
% error('debugging in bry_subgrid')
 end
end

