
      chi_grd = '/batavia/WEDDY/weddy_grd.nc';
      par_grd = '/batavia/GULFS/gulfs_grd.nc';


      nc    = netcdf(chi_grd, 'nowrite');
      lonc  = nc{'lon_rho'}(:);
      latc  = nc{'lat_rho'}(:);
      close(nc)

      nc    = netcdf(par_grd, 'nowrite');
      lonp  = nc{'lon_rho'}(:);
      latp  = nc{'lat_rho'}(:);
      [Mpp,Lpp] = size(latp);
      close(nc)

      display('going delaunay');
      tri_fullpar = delaunay(lonp,latp);
      display('returndelaunay');

    % Compute minimal subgrid extracted from full parent grid
      t = squeeze(tsearch(lonp,latp,tri_fullpar,lonc,latc));
    % Deal with child points that are outside parent grid (those points should be masked!)
      if (length(t(~isfinite(t)))>0);
       disp('Warning in new_bry_subgrid: outside point(s) detected.');
       [lonc,latc] = fix_outside_child(lonc,latc,t);
       t  = squeeze(tsearch(lonp,latp,tri_fullpar,lonc,latc));
      end;
      index       = tri_fullpar(t,:);
      [idxj,idxi] = ind2sub([Mpp Lpp], index);

      imin = min(min(idxi));
      imax = max(max(idxi));
      jmin = min(min(idxj));
      jmax = max(max(idxj));

      [imin imax jmin jmax]

      plot(lonp(jmin:jmax,imin:imax),latp(jmin:jmax,imin:imax),'.k')
      hold on
      plot(lonc,latc,'.r')
      hold off

