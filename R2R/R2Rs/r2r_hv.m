function roms2roms_hv(par_grd,par_data,par_tind, chd_grd,       ...
                     bry_name, chd_tind, chdscd,   parscd,      ...
                     obcflag, limits)
%--------------------------------------------------------------
%  Extract arbitrary variables from parent-ROMS history file
%  to save in boundary perimeter files on child-ROMS grid.
%
%  input
%  =====
%  par_grd:         parent grid .nc file name
%  par_data:        parent history .nc file name
%  par_tind:        time index for parent history file
%  chd_grd:         child grid .nc file name
%  bry_name:        child boundary .nc file name
%  chd_tind:        time index for child boundary file
%  chdscd:          child grid s-coordinate parameters (object)
%  parscd:          parent grid s-coordinate parameters (object)
%  obcflag:         open boundaries flag (1=open , [S E N W])
%
%  Inspired by Roms_tools (IRD).
%  Thanks to Pierrick, Patrick (IRD), Francois (UCLA), Yusuke (UCLA)
%  Heavily modified from French produce
%  Jeroen Molemaker (UCLA); nmolem@ucla.edu
%    Evan Mason (ULPGC), June 2007
%--------------------------------------------------------------
%

% Get S-coordinate params for child grid
  theta_b_c = chdscd.theta_b;
  theta_s_c = chdscd.theta_s;
  hc_c      = chdscd.hc;
  N_c       = chdscd.N;
  scoord_c  = chdscd.scoord;

% Get S-coordinate params for parent grid
  theta_b_p = parscd.theta_b;
  theta_s_p = parscd.theta_s;
  hc_p      = parscd.hc;
  N_p       = parscd.N;
  scoord_p  = parscd.scoord;

  nc = netcdf(chd_grd, 'nowrite');
  [mpc npc] = size(nc{'h'}(:));
  close(nc)

  nd   = netcdf(bry_name, 'write');

% Set bry_time
  np = netcdf(par_data, 'nowrite');
  time = np{'ocean_time'}(par_tind);
  par_data
  np{'ocean_time'}(:) 
  pind = 1;
  for tind = chd_tind
     [tind pind]
    nd{'bry_time'}(tind) = time(pind)/(3600*24.);
    pind = pind + 1;
  end
  close(np)

 for bnd = 1:4
  disp('-------------------------------------------------------------')
  if ~obcflag(bnd)
    disp('Closed boundary')
    continue
  end
  if bnd==1 
   disp('South boundary')
   i0 =  1; 
   i1 = npc;
   j0 = 1;
   j1 = 2;
   fcoef = 'r2r_coefs_south.mat';
  end
  if bnd==2 
   disp('East boundary')
   i0 = npc-1;
   i1 = npc;
   j0 = 1;
   j1 = mpc;
   fcoef = 'r2r_coefs_east.mat';
  end
  if bnd==3 
   disp('North boundary')
   i0 = 1;
   i1 = npc;
   j0 = mpc-1;
   j1 = mpc;
   fcoef = 'r2r_coefs_north.mat';
  end
  if bnd==4 
   disp('West boundary')
   i0 = 1;
   i1 = 2;
   j0 = 1;
   j1 = mpc;
   fcoef = 'r2r_coefs_west.mat';
  end

% Compute minimal subgrid extracted from parent grid
  imin = limits(bnd,1); imax = limits(bnd,2);
  jmin = limits(bnd,3); jmax = limits(bnd,4);

% Get topography data from childgrid
  nc = netcdf(chd_grd, 'nowrite');
  hc    = nc{'h'}(j0:j1,i0:i1);
  pm    = nc{'pm'}(j0:j1,i0:i1);       pn    = nc{'pn'}(j0:j1,i0:i1);
  maskc = nc{'mask_rho'}(j0:j1,i0:i1); angc  = nc{'angle'}(j0:j1,i0:i1);
  lonc  = nc{'lon_rho'}(j0:j1,i0:i1);  latc  = nc{'lat_rho'}(j0:j1,i0:i1);
  cosc  = cos(angc);         sinc  = sin(angc);

  tra_lon = nc{'tra_lon'}(:);
  tra_lat = nc{'tra_lat'}(:);
  rotate  = nc{'rotate'}(:);
  xy_flip = nc{'xy_flip'}(:);

  close(nc);


  [Mc,Lc] = size(maskc);
  maskc3d = zeros(N_c,Mc,Lc);
  for k = 1:N_c
   maskc3d(k,:,:) = maskc;
  end
  umask = maskc3d(:,:,2:end).*maskc3d(:,:,1:end-1);
  vmask = maskc3d(:,2:end,:).*maskc3d(:,1:end-1,:);

% Get parent grid and squeeze minimal subgrid
  np    = netcdf(par_grd, 'nowrite');
  masks = squeeze(np{'mask_rho'}(jmin:jmax, imin:imax));
  lons  = squeeze(np{'lon_rho'}(jmin:jmax, imin:imax));%  + 360; 
  lats  = squeeze(np{'lat_rho'}(jmin:jmax, imin:imax));
  angs  = squeeze(np{'angle'}(jmin:jmax, imin:imax));
  hs    = squeeze(np{'h'}(jmin:jmax, imin:imax));
  close(np);
  coss = cos(angs); sins = sin(angs);
  if sum(isnan(masks))>0
      disp('Setting NaNs in masks to zero')
      masks(isnan(masks))=0;
      disp('You probably have land masking defined in cppdefs.h...')
  end

  % Z-coordinate (3D) on minimal subgrid and child grid
  % Sasha recommends the 0 multiplication
  zs = zlevs3(hs, hs*0, theta_s_p, theta_b_p, hc_p, N_p, 'r', scoord_p);
  zc = zlevs3(hc, hc*0, theta_s_c, theta_b_c, hc_c, N_c, 'r', scoord_c);
  zw = zlevs3(hc, hc*0, theta_s_c, theta_b_c, hc_c, N_c, 'w', scoord_c);

  [Np Mp Lp] = size(zs);
  [Nc Mc Lc] = size(zc);

  if exist(fcoef,'file')
    disp(' ')
    disp('Reading interpolation coefficients from file');
    load(fcoef)
  else
    tic
    disp(' ')
    disp('Computing interpolation coefficients');
    [elem2d,coef2d,nnel] = get_tri_coef(lons,lats,lonc,latc,masks);

    A = get_hv_coef(zs, zc, coef2d, elem2d, lons, lats, lonc, latc);
    save(fcoef,'elem2d','coef2d','nnel','A')
    toc
  end

% Open parent data file
  np = netcdf(par_data, 'nowrite');
  [ntp,nyp,nxp] = size(np{'zeta'}(:));

% Loop on parent time index
  loop = 1;
  for tind = par_tind
    tout = chd_tind(loop);
    loop = loop + 1;
%   disp('-------------------------------------------------------------')
%   disp(['tind = ' int2str(tind) ' of ' int2str(length(par_tind))...
%         ' steps. (tout=' int2str(tout) ')']);

%   Surface elevation on minimal subgrid and child grid
%   disp('--- zeta')
    zetas = squeeze(np{'zeta'}(tind,jmin:jmax,imin:imax));
    zetas = fillmask(zetas, 1, masks, nnel);
    zetac = sum(coef2d.*zetas(elem2d), 3);

%   Prepare for estimating barotropic velocity
    dz  = zw(2:end,:,:)-zw(1:end-1,:,:);
    dzu = 0.5*(dz(:,:,1:end-1)+dz(:,:,2:end));
    dzv = 0.5*(dz(:,1:end-1,:)+dz(:,2:end,:));

%   Process scalar 3D variables
    for vint = 1:2 % Loop on the tracers
      if (vint==1)
        svar='temp';
      elseif (vint==2)
        svar='salt';
      end
%     disp(['--- ' svar])
      var = squeeze(np{svar}(tind,:,jmin:jmax,imin:imax));
      var = fillmask(var,1,masks,nnel);
      var = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);
      if (bnd ==1 )
        nd{[svar '_south']}(tout,:,:)= var(:, 1, :);
      end
      if (bnd ==2 )
        nd{[svar '_east']}(tout,:,:) = var(:, :,end);
      end
      if (bnd == 3 )
        nd{[svar '_north']}(tout,:,:)= var(:,end, :);
      end
      if (bnd == 4 )
        nd{[svar '_west']}(tout,:,:) = var(:, :, 1);
      end
    end  % End loop on vint

    %  Read in staggered velocities and move to rho-points
    if (imin==1|imax==nxp)
     ud = squeeze(np{'u'}(tind,:,jmin:jmax  ,imin:imax-1));
     ur = u2rho(ud);
    else
     ud = squeeze(np{'u'}(tind,:,jmin:jmax  ,imin-1:imax));
     ur = 0.5*(ud(:,:,1:end-1) + ud(:,:,2:end));
    end
    if (jmin==1|jmax==nyp)
     vd = squeeze(np{'v'}(tind,:,jmin:jmax-1,imin  :imax));
     vr = v2rho(vd);
    else
     vd = squeeze(np{'v'}(tind,:,jmin-1:jmax,imin  :imax));
     vr = 0.5*(vd(:,1:end-1,:) + vd(:,2:end,:));
    end

    % Rotate to north
    us = zeros(Np,Mp,Lp);
    vs = zeros(Np,Mp,Lp);
    for k = 1:Np
      us(k,:,:) = squeeze(ur(k,:,:)).*coss - squeeze(vr(k,:,:)).*sins;
      vs(k,:,:) = squeeze(vr(k,:,:)).*coss + squeeze(ur(k,:,:)).*sins;
    end

    % 3d interpolation of us and vs to child grid

    us = fillmask(us, 0, masks, nnel);
    vs = fillmask(vs, 0, masks, nnel);
    ud = reshape(A*reshape(us, Np*Mp*Lp,1), Nc,Mc,Lc);
    vd = reshape(A*reshape(vs, Np*Mp*Lp,1), Nc,Mc,Lc);

    % Rotate to child orientation
    us = zeros(Nc, Mc, Lc);
    vs = zeros(Nc, Mc, Lc);
    for k=1:Nc
      us(k,:,:) = squeeze(ud(k,:,:)).*cosc + squeeze(vd(k,:,:)).*sinc;
      vs(k,:,:) = squeeze(vd(k,:,:)).*cosc - squeeze(ud(k,:,:)).*sinc;
    end

    %% back to staggered u and v points
    u = 0.5*(us(:,:,1:Lc-1) + us(:,:,2:Lc));
    v = 0.5*(vs(:,1:Mc-1,:) + vs(:,2:Mc,:));

    u = u.*umask;
    v = v.*vmask;

    % Get barotropic velocity
    if sum(isnan(u)) > 0
      error('nans in u velocity!')
    end
    if sum(isnan(v)) > 0
      error('nans in v velocity!')
    end

%   disp('--- velocities');
    hu   = sum(dzu.*u); hv   = sum(dzv.*v);
    D_u  = sum(dzu);    D_v  = sum(dzv);
    [dum Mu Lu] = size(hu);
    [dum Mv Lv] = size(hv);
    ubar = reshape(hu./D_u,Mu, Lu);     
    vbar = reshape(hv./D_v,Mv, Lv);


    % Save perimeter zeta, ubar, vbar, u and v data to bryfile
%   disp(['------ saving velocity fields for record ' ...
%                       int2str(tout) ' to bry file'])
    if (bnd ==1)
      nd{'ubar_south'}(tout,:)= squeeze(ubar(1,:));
      nd{'vbar_south'}(tout,:)= squeeze(vbar(1,:));
      nd{'zeta_south'}(tout,:)= squeeze(zetac(1,:));
      nd{'u_south'}(tout,:,:) = squeeze(u(:,1,:));
      nd{'v_south'}(tout,:,:) = squeeze(v(:,1,:));
    end
    if (bnd == 2)
      nd{'ubar_east'}(tout,:) = squeeze(ubar(:,Lu));
      nd{'vbar_east'}(tout,:) = squeeze(vbar(:,Lv));
      nd{'zeta_east'}(tout,:) = squeeze(zetac(:,Lc));
      nd{'u_east'}(tout,:,:)  = squeeze(u(:,:,Lu));
      nd{'v_east'}(tout,:,:)  = squeeze(v(:,:,Lv));
    end
    if (bnd == 3)
      nd{'ubar_north'}(tout,:)= squeeze(ubar(Mu,:));
      nd{'vbar_north'}(tout,:)= squeeze(vbar(Mv,:));
      nd{'zeta_north'}(tout,:)= squeeze(zetac(Mc,:));
      nd{'u_north'}(tout,:,:) = squeeze(u(:,Mu,:));
      nd{'v_north'}(tout,:,:) = squeeze(v(:,Mv,:));
    end
    if (bnd == 4)
      nd{'ubar_west'}(tout,:) = squeeze(ubar(:,1));
      nd{'vbar_west'}(tout,:) = squeeze(vbar(:,1));
      nd{'zeta_west'}(tout,:) = squeeze(zetac(:,1));
      nd{'u_west'}(tout,:,:)  = squeeze(u(:,:,1));
      nd{'v_west'}(tout,:,:)  = squeeze(v(:,:,1));
    end

    % Close child bryfile

  end    % End loop on par_tind

 end    % End loop bnd

 close(nd)
 close(np)

  return























