function limits = r2r_subgrid(parentgrid,childgrid,obcflag);
%
%   Find lower and upper index in i and j for the minimal
%   parent grid that contains all of the boundary bnd of
%   the child grid.
%
%   (c) 2007,  Jeroen Molemaker
%
%--------------------------------------------------------------

% get data from parentgrid
 nc=netcdf(parentgrid,'nowrite');
 lonp = nc{'lon_rho'}(:);%  +360;
 latp = nc{'lat_rho'}(:);
 tra_lon = nc{'tra_lon'}(:);
 tra_lat = nc{'tra_lat'}(:);
 rotate  = nc{'rotate'}(:);
 xy_flip = nc{'xy_flip'}(:);
 close(nc);

 [lonp,latp] = grd2mesh(lonp,latp,tra_lon,tra_lat,rotate,xy_flip,1);
 [ny,nx] = size(lonp);

 limits = zeros(4,4);
 for bnd = 1:4
  if ~obcflag(bnd)
     continue
  end
  nc = netcdf(childgrid,'nowrite');
  if bnd == 1 %% South
    lonc = nc{'lon_rho'}(1:2,1:end);
    latc = nc{'lat_rho'}(1:2,1:end);
  end
  if bnd == 2 %% East
    lonc = nc{'lon_rho'}(1:end,end-1:end); 
    latc = nc{'lat_rho'}(1:end,end-1:end);
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
 
 
   [lonc,latc] = grd2mesh(lonc,latc,tra_lon,tra_lat,rotate,xy_flip,0);

   if min(min(lonc))<lonp(1,1)
    i0 = 1;
   else
    i0 = find(lonp(1,:)<min(min(lonc)),1,'last');
   end
   if max(max(lonc))>lonp(1,nx)
    i0 = nx;
   else
    i1 = find(lonp(1,:)>max(max(lonc)),1,'first');
   end
   if min(min(latc))<latp(1,1)
    j0 = 1;
   else
    j0 = find(latp(:,1)<min(min(latc)),1,'last');
   end
   if max(max(latc))>latp(ny,1)
     j1 = ny;
   else
     j1 = find(latp(:,1)>max(max(latc)),1,'first') ;
   end


   limits(bnd,1)  = i0;
   limits(bnd,2)  = i1;
   limits(bnd,3)  = j0;
   limits(bnd,4)  = j1;
 
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
    limits(bnd,:);
%  error('debugging in bry_subgrid')
  end
end

