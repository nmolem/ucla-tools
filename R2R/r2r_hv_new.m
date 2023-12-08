function roms2roms_hv(par_grd,par_data,par_tind, chd_grd,       ...
                     bry_name, chd_tind, chdscd,   parscd,      ...
                     obcflag, limits)
  persistent elem2d_s coef2d_s nnel_s A_s elem2d_e coef2d_e nnel_e A_e  ...
             elem2d_n coef2d_n nnel_n A_n elem2d_w coef2d_w nnel_w A_w
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
  par_file  = parscd.par_file;

% nc = netcdf(chd_grd, 'nowrite');
% [mpc npc] = size(nc{'h'}(:));
% close(nc)
  [mpc npc] = size(ncread(chd_grd,'h')');

% nd   = netcdf(bry_name, 'write');

% Set bry_time
% np = netcdf(par_data, 'nowrite');
% time = np{'ocean_time'}(par_tind);
% np{'ocean_time'}(:) 
  time = ncread(par_data,'ocean_time',par_tind,1);
  pind = 1;
  for tind = chd_tind
     [tind pind]
%   nd{'bry_time'}(tind) = time(pind)/(3600*24.);
    ncwrite(bry_name,'bry_time',time(pind)/(3600*24.),tind);
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
% nc = netcdf(chd_grd, 'nowrite');
% hc    = nc{'h'}(j0:j1,i0:i1);
% pm    = nc{'pm'}(j0:j1,i0:i1);       pn    = nc{'pn'}(j0:j1,i0:i1);
% maskc = nc{'mask_rho'}(j0:j1,i0:i1); angc  = nc{'angle'}(j0:j1,i0:i1);
% lonc  = nc{'lon_rho'}(j0:j1,i0:i1);  latc  = nc{'lat_rho'}(j0:j1,i0:i1);
% close(nc);

  h    = ncread(chd_grd,'h',[i0 j0],[i1-i0+1 j1-j0+1])';
  pm   = ncread(chd_grd,'pm',[i0 j0],[i1-i0+1 j1-j0+1])';
  pn   = ncread(chd_grd,'pn',[i0 j0],[i1-i0+1 j1-j0+1])';
  maskc= ncread(chd_grd,'mask_rho',[i0 j0],[i1-i0+1 j1-j0+1])';
  angc = ncread(chd_grd,'angle',[i0 j0],[i1-i0+1 j1-j0+1])';
  lonc = ncread(chd_grd,'lon_rho',[i0 j0],[i1-i0+1 j1-j0+1])';
  latc = ncread(chd_grd,'lat_rho',[i0 j0],[i1-i0+1 j1-j0+1])';

  cosc  = cos(angc);         sinc  = sin(angc);

  [Mc,Lc] = size(maskc);
  maskc3d = zeros(N_c,Mc,Lc);
  for k = 1:N_c
   maskc3d(k,:,:) = maskc;
  end
  umask = maskc3d(:,:,2:end).*maskc3d(:,:,1:end-1);
  vmask = maskc3d(:,2:end,:).*maskc3d(:,1:end-1,:);

% Get parent grid and squeeze minimal subgrid
% np    = netcdf(par_grd, 'nowrite');
% masks = squeeze(np{'mask_rho'}(jmin:jmax, imin:imax));
% lons  = squeeze(np{'lon_rho'}(jmin:jmax, imin:imax)); 
% lats  = squeeze(np{'lat_rho'}(jmin:jmax, imin:imax));
% angs  = squeeze(np{'angle'}(jmin:jmax, imin:imax));
% hs    = squeeze(np{'h'}(jmin:jmax, imin:imax));
% close(np);

  masks= ncread(chd_grd,'mask_rho',[imin jmin],[imax-imin+1 jmax-jmin+1])';
  angs = ncread(chd_grd,'angle',[imin jmin],[imax-imin+1 jmax-jmin+1])';
  lons = ncread(chd_grd,'lon_rho',[imin jmin],[imax-imin+1 jmax-jmin+1])';
  lats = ncread(chd_grd,'lat_rho',[imin jmin],[imax-imin+1 jmax-jmin+1])';
  hs   = ncread(chd_grd,'h',[imin jmin],[imax-imin+1 jmax-jmin+1])';

  coss = cos(angs); sins = sin(angs);
  if sum(isnan(masks))>0
      disp('Setting NaNs in masks to zero')
      masks(isnan(masks))=0;
      disp('You probably have land masking defined in cppdefs.h...')
  end

  % Z-coordinate (3D) on minimal subgrid and child grid
  % Sasha recommends the 0 multiplication
% zs = zlevs3(hs, hs*0, theta_s_p, theta_b_p, hc_p, N_p, 'r', scoord_p);
  zs = zlev_cf(hs, hs*0, par_grd, 'r');
  zc = zlevs3(hc, hc*0, theta_s_c, theta_b_c, hc_c, N_c, 'r', scoord_c);
  zw = zlevs3(hc, hc*0, theta_s_c, theta_b_c, hc_c, N_c, 'w', scoord_c);

  [Np Mp Lp] = size(zs);
  [Nc Mc Lc] = size(zc);

  if exist(fcoef,'file')
    disp(' ')
    disp('Reading interpolation coefficients from file');
%   load(fcoef)
    if bnd==1
      elem2d = elem2d_s;
      coef2d = coef2d_s;
      nnel   = nnel_s;
      A      = A_s;
    end
    if bnd==2
      elem2d = elem2d_e;
      coef2d = coef2d_e;
      nnel   = nnel_e;
      A      = A_e;
    end
    if bnd==3
      elem2d = elem2d_n;
      coef2d = coef2d_n;
      nnel   = nnel_n;
      A      = A_n;
    end
    if bnd==4
      elem2d = elem2d_w;
      coef2d = coef2d_w;
      nnel   = nnel_w;
      A      = A_w;
    end
  else
    tic
    disp(' ')
    disp('Computing interpolation coefficients');
    [elem2d,coef2d,nnel] = get_tri_coef(lons,lats,lonc,latc,masks);

    A = get_hv_coef(zs, zc, coef2d, elem2d, lons, lats, lonc, latc);
    if bnd==1
      elem2d_s = elem2d;
      coef2d_s = coef2d;
      nnel_s   = nnel;
      A_s      = A;
    end
    if bnd==2
      elem2d_e = elem2d;
      coef2d_e = coef2d;
      nnel_e   = nnel;
      A_e      = A;
    end
    if bnd==3
      elem2d_n = elem2d;
      coef2d_n = coef2d;
      nnel_n   = nnel;
      A_n      = A;
    end
    if bnd==4
      elem2d_w = elem2d;
      coef2d_w = coef2d;
      nnel_w   = nnel;
      A_w      = A;
    end
    save(fcoef,'elem2d','coef2d','nnel','A')
    toc
  end

% Open parent data file
% np = netcdf(par_data, 'nowrite');
% [ntp,nyp,nxp] = size(np{'zeta'}(:));
  [ntp,nyp,nxp] = size(ncread(par_data,'zeta'));

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
%     size(var)
%     [Np Mp Lp]
      var = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);
%     size(var)
%     size(nd{[svar '_south']}(tout,:,:))
%     disp('afmetingen')
      if (bnd ==1 )
%       nd{[svar '_south']}(tout,:,:)= var(:, 1, :);
	ncwrite(bry_name,[svar '_south'],squeeze(var(:,1,:))',[:,:,tout]);
      end
      if (bnd ==2 )
%       nd{[svar '_east']}(tout,:,:) = var(:, :,end);
	ncwrite(bry_name,[svar '_east'],squeeze(var(:,:,end))',[:,:,tout]);
      end
      if (bnd == 3 )
%       nd{[svar '_north']}(tout,:,:)= var(:,end, :);
	ncwrite(bry_name,[svar '_north'],squeeze(var(:,end,:))',[:,:,tout]);
      end
      if (bnd == 4 )
%       nd{[svar '_west']}(tout,:,:) = var(:, :, 1);
	ncwrite(bry_name,[svar '_west'],squeeze(var(:,:,1))',[:,:,tout]);
      end
    end  % End loop on vint

    %  Read in staggered velocities and move to rho-points
    if (imin==1|imax==nxp)
%    ud = squeeze(np{'u'}(tind,:,jmin:jmax  ,imin:imax-1));ud(abs(ud)>1e2) = 0;
     ud = squeeze(ncread(par_data,'u',[imin jmin 1 tind],[imax-imin jmax-jmin+1 inf 1]));
     ud = permute(ud,[3 2 1]);
     ud(abs(ud)>1e2) = 0;
     ur = u2rho(ud);
    else
%    ud = squeeze(np{'u'}(tind,:,jmin:jmax  ,imin-1:imax));ud(abs(ud)>1e2) = 0;
     ud = squeeze(ncread(par_data,'u',[imin-1 jmin 1 tind],[imax-imin+2 jmax-jmin+1 inf 1]));
     ud = permute(ud,[3 2 1]);
     ud(abs(ud)>1e2) = 0;
     ur = 0.5*(ud(:,:,1:end-1) + ud(:,:,2:end));
    end
    if (jmin==1|jmax==nyp)
%    vd = squeeze(np{'v'}(tind,:,jmin:jmax-1,imin  :imax));vd(abs(vd)>1e2) = 0;
     vd = squeeze(ncread(par_data,'v',[imin jmin 1 tind],[imax-imin+1 jmax-jmin inf 1]));
     vd = permute(vd,[3 2 1]);
     vd(abs(vd)>1e2) = 0;
     vr = v2rho(vd);
    else
%    vd = squeeze(np{'v'}(tind,:,jmin-1:jmax,imin  :imax));vd(abs(vd)>1e2) = 0;
     vd = squeeze(ncread(par_data,'v',[imin jmin-1 1 tind],[imax-imin+1 jmax-jmin+2 inf 1]));
     vd = permute(vd,[3 2 1]);
     vd(abs(vd)>1e2) = 0;
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
%     nd{'ubar_south'}(tout,:)= squeeze(ubar(1,:));
%     nd{'vbar_south'}(tout,:)= squeeze(vbar(1,:));
%     nd{'zeta_south'}(tout,:)= squeeze(zetac(1,:));
%     nd{'u_south'}(tout,:,:) = squeeze(u(:,1,:));
%     nd{'v_south'}(tout,:,:) = squeeze(v(:,1,:));
      ncwrite(bry_name,'ubar_south',squeeze(ubar(1,:))',[1 tout]);
      ncwrite(bry_name,'vbar_south',squeeze(vbar(1,:))',[1 tout]);
      ncwrite(bry_name,'zeta_south',squeeze(zetac(1,:))',[1 tout]);
      ncwrite(bry_name,'u_south'   ,squeeze(u(:,1,:))',[1 1 tout]);
      ncwrite(bry_name,'v_south'   ,squeeze(v(:,1,:))',[1 1 tout]);
    end
    if (bnd == 2)
%     nd{'ubar_east'}(tout,:) = squeeze(ubar(:,Lu));
%     nd{'vbar_east'}(tout,:) = squeeze(vbar(:,Lv));
%     nd{'zeta_east'}(tout,:) = squeeze(zetac(:,Lc));
%     nd{'u_east'}(tout,:,:)  = squeeze(u(:,:,Lu));
%     nd{'v_east'}(tout,:,:)  = squeeze(v(:,:,Lv));
      ncwrite(bry_name,'ubar_east',squeeze(ubar(:,Lu))',[1 tout]);
      ncwrite(bry_name,'vbar_east',squeeze(vbar(:,Lv))',[1 tout]);
      ncwrite(bry_name,'zeta_east',squeeze(zetac(:,Lc))',[1 tout]);
      ncwrite(bry_name,'u_east'   ,squeeze(u(:,:,Lu))',[1 1 tout]);
      ncwrite(bry_name,'v_east'   ,squeeze(v(:,:,Lv))',[1 1 tout]);
    end
    if (bnd == 3)
%     nd{'ubar_north'}(tout,:)= squeeze(ubar(Mu,:));
%     nd{'vbar_north'}(tout,:)= squeeze(vbar(Mv,:));
%     nd{'zeta_north'}(tout,:)= squeeze(zetac(Mc,:));
%     nd{'u_north'}(tout,:,:) = squeeze(u(:,Mu,:));
%     nd{'v_north'}(tout,:,:) = squeeze(v(:,Mv,:));
      ncwrite(bry_name,'ubar_north',squeeze(ubar(Mu,:))',[1 tout]);
      ncwrite(bry_name,'vbar_north',squeeze(vbar(Mv,:))',[1 tout]);
      ncwrite(bry_name,'zeta_north',squeeze(zetac(Mc,:))',[1 tout]);
      ncwrite(bry_name,'u_north'   ,squeeze(u(:,Mu,:))',[1 1 tout]);
      ncwrite(bry_name,'v_north'   ,squeeze(v(:,Mv,:))',[1 1 tout]);
    end
    if (bnd == 4)
%     nd{'ubar_west'}(tout,:) = squeeze(ubar(:,1));
%     nd{'vbar_west'}(tout,:) = squeeze(vbar(:,1));
%     nd{'zeta_west'}(tout,:) = squeeze(zetac(:,1));
%     nd{'u_west'}(tout,:,:)  = squeeze(u(:,:,1));
%     nd{'v_west'}(tout,:,:)  = squeeze(v(:,:,1));
      ncwrite(bry_name,'ubar_west',squeeze(ubar(:,1))',[1 tout]);
      ncwrite(bry_name,'vbar_west',squeeze(vbar(:,1))',[1 tout]);
      ncwrite(bry_name,'zeta_west',squeeze(zetac(:,1))',[1 tout]);
      ncwrite(bry_name,'u_west'   ,squeeze(u(:,:,1))',[1 1 tout]);
      ncwrite(bry_name,'v_west'   ,squeeze(v(:,:,1))',[1 1 tout]);
    end

    % Close child bryfile

  end    % End loop on par_tind

 end    % End loop bnd

 close(nd)
 close(np)
 close all
 fclose('all')

  return























