function r2r_hv(par_grd,par_data,par_tind, chd_grd,       ...
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

  [mpc npc] = size(ncread(chd_grd,'h')');

% Set bry_time
% np = netcdf(par_data, 'nowrite');
% np{'ocean_time'}(:) 
% time = np{'ocean_time'}(par_tind);
% disp('r2r_hv')
% par_tind
% time = ncread(par_data,'ocean_time',[par_tind],1);
  time = ncread(par_data,'ocean_time');
  pind = 1;
  for tind = chd_tind
%   [tind pind]
    ncwrite(bry_name,'bry_time',time(pind)/(3600*24.),[tind]);
%   nd{'bry_time'}(tind) = time(pind)/(3600*24.);
%   nd{'bry_time'}(tind) = time(pind)/(3600*24.) - 55197;
    pind = pind + 1;
  end

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
  hc = ncread(chd_grd,'h')';  hc = hc(j0:j1,i0:i1);
  pm = ncread(chd_grd,'pm')'; pm = pm(j0:j1,i0:i1);
  pn = ncread(chd_grd,'pn')'; pn = pn(j0:j1,i0:i1);
  maskc = ncread(chd_grd,'mask_rho')'; maskc = maskc(j0:j1,i0:i1);
  angc  = ncread(chd_grd,'angle')';    angc  = angc(j0:j1,i0:i1);
  lonc  = ncread(chd_grd,'lon_rho')';  lonc  = lonc(j0:j1,i0:i1);
  latc  = ncread(chd_grd,'lat_rho')';  latc  = latc(j0:j1,i0:i1);
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
 [jmin jmax imin imax]
  size(ncread(par_grd,'mask_rho'))

  masks = ncread(par_grd,'mask_rho')'; masks = masks(jmin:jmax,imin:imax);
  lons  = ncread(par_grd,'lon_rho' )'; lons  =  lons(jmin:jmax,imin:imax);
  lats  = ncread(par_grd,'lat_rho' )'; lats  =  lats(jmin:jmax,imin:imax);
  angs  = ncread(par_grd,'angle'   )'; angs  =  angs(jmin:jmax,imin:imax);
  hs    = ncread(par_grd,'h'       )'; hs    =    hs(jmin:jmax,imin:imax);
  coss = cos(angs); sins = sin(angs);
  lons  = double(lons);
  lats  = double(lats);

  if sum(isnan(masks))>0
      disp('Setting NaNs in masks to zero')
      masks(isnan(masks))=0;
      disp('You probably have land masking defined in cppdefs.h...')
  end

  % Z-coordinate (3D) on minimal subgrid and child grid
  % Sasha recommends the 0 multiplication
  zs = zlevs3(hs, hs*0, theta_s_p, theta_b_p, hc_p, N_p, 'r', scoord_p);
% zs = zlev_cf(hs, hs*0, par_grd, 'r');
  zc = zlevs3(hc, hc*0, theta_s_c, theta_b_c, hc_c, N_c, 'r', scoord_c);
  zw = zlevs3(hc, hc*0, theta_s_c, theta_b_c, hc_c, N_c, 'w', scoord_c);

  [Np Mp Lp] = size(zs);
  [Nc Mc Lc] = size(zc);

  if exist(fcoef,'file')
    disp(' ')
    disp('Reading interpolation coefficients from file');
    load(fcoef)
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

  [nxp,nyp,ntp] = size(ncread(par_data,'zeta'));

% Loop on parent time index
  loop = 1;
  for tind = par_tind
    tout = chd_tind(loop);
    loop = loop + 1;
    disp('-------------------------------------------------------------')
    disp(['tind = ' int2str(tind) ' of ' int2str(length(par_tind))...
          ' steps. (tout=' int2str(tout) ')']);

%   Surface elevation on minimal subgrid and child grid

%   disp('--- zeta')
    zetas = ncread(par_data,'zeta'  ,[imin jmin tind],[imax-imin+1,jmax-jmin+1 1])';
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
      var = squeeze(ncread(par_data,svar,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = double(permute(var,[3 2 1]));
      var = fillmask(var,1,masks,nnel);
      var = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      if (bnd ==1 )
        ncwrite(bry_name,[svar '_south'],squeeze(var(:,  1,:))',[1 1 tout]);
      end
      if (bnd ==2 )
        ncwrite(bry_name,[svar '_east'] ,squeeze(var(:,:,end))',[1 1 tout]);
      end
      if (bnd == 3 )
        ncwrite(bry_name,[svar '_north'],squeeze(var(:,end,:))',[1 1 tout]);
      end
      if (bnd == 4 )
        ncwrite(bry_name,[svar '_west'] ,squeeze(var(:,:,  1))',[1 1 tout]);
      end
    end  % End loop on vint

    %  Read in staggered velocities and move to rho-points
    if (imin==1|imax==nxp)
     ud = squeeze(ncread(par_data,'u',[imin jmin 1 tind],[imax-imin,jmax-jmin+1 inf 1]));
     ud = double(permute(ud,[3 2 1]));
     ud(abs(ud)>1e2) = 0;
     ur = u2rho(ud);
    else
     ud = squeeze(ncread(par_data,'u',[imin-1 jmin 1 tind],[imax-imin+2,jmax-jmin+1 inf 1]));
     ud = double(permute(ud,[3 2 1]));
     ud(abs(ud)>1e2) = 0;
     ur = 0.5*(ud(:,:,1:end-1) + ud(:,:,2:end));
    end
    if (jmin==1|jmax==nyp)
     vd = squeeze(ncread(par_data,'v',[imin jmin 1 tind],[imax-imin+1,jmax-jmin inf 1]));
     vd = double(permute(vd,[3 2 1]));
     vd(abs(vd)>1e2) = 0;
     vr = v2rho(vd);
    else
     vd = squeeze(ncread(par_data,'v',[imin jmin-1 1 tind],[imax-imin+1,jmax-jmin+2 inf 1]));
     vd = double(permute(vd,[3 2 1]));
     vd(abs(vd)>1e2) = 0;
     vr = 0.5*(vd(:,1:end-1,:) + vd(:,2:end,:));
    end

    % consider cases like jmin==1 -> add an extra row of w-points
    % Also, the code above seems odd, jmin==1 and jmax==nyp should be considered separately
    % so, for now, just flag these cases for a warning
    if (jmin==1|jmax==nyp|imin==1|imax==nxp)
       display('jmin==1|jmax==nyp|imin==1|imax==nxp not properly implemented yet')
%      error
    end
    % compute vertical flux on parent grid
%   hse = ncread(par_grd,'h',[imin-1 jmin-1],[imax-imin+3,jmax-jmin+3])';
%   pms = ncread(par_grd,'pm',[imin-1 jmin-1],[imax-imin+3,jmax-jmin+3])';
%   pns = ncread(par_grd,'pn',[imin-1 jmin-1],[imax-imin+3,jmax-jmin+3])';
%   zws = zlevs3(hse, hse*0, theta_s_p, theta_b_p, hc_p, N_p, 'w', scoord_p);
%   dzs = zws(2:end,:,:)-zws(1:end-1,:,:);
%   dzus = 0.5*(dzs(:,2:end-1,1:end-1)+dzs(:,2:end-1,2:end));
%   dzvs = 0.5*(dzs(:,1:end-1,2:end-1)+dzs(:,2:end,2:end-1));

%   pnus = 0.5*(pns(2:end-1,1:end-1)+pns(2:end-1,2:end));
%   pmvs = 0.5*(pms(1:end-1,2:end-1)+pms(2:end,2:end-1));

%   size(pnus)
%   size(dzus)
%   size(ud)
%   Flxu = ud;
%   Flxv = vd;

%   for k = 1:Np
%     Flxu(k,:,:) = squeeze(dzus(k,:,:).*ud(k,:,:))./pnus;
%     Flxv(k,:,:) = squeeze(dzvs(k,:,:).*vd(k,:,:))./pmvs;
%   end

    ws = zeros(Np+1,Mp,Lp);
%   for k = 1:Np
%    ws(k+1,:,:) = squeeze(ws(k,:,:)) - pms(2:end-1,2:end-1).*pns(2:end-1,2:end-1).*squeeze( ...
%               Flxu(k,:,2:end)-Flxu(k,:,1:end-1) +  ...
%               Flxv(k,2:end,:)-Flxv(k,1:end-1,:) );
%   end
    % average ws to rho points
    ws = 0.5*(ws(2:end,:,:)+ws(1:end-1,:,:));


    if 0 % compute vertical velocity
      zse = zlevs3(hse, hse*0, theta_s_p, theta_b_p, hc_p, N_p, 'r', scoord_p);
      dzdx = 0.5*(zse(:,2:end-1,2:end)-zse(:,2:end-1,1:end-1));
      dzdy = 0.5*(zse(:,2:end,2:end-1)-zse(:,1:end-1,2:end-1));
      ws = ws + 0.5*( dzdx(:,:,2:end).*ud(:,:,2:end) + dzdx(:,:,1:end-1).*ud(:,:,1:end-1) ) ...
              + 0.5*( dzdy(:,2:end,:).*vd(:,2:end,:) + dzdy(:,1:end-1,:).*vd(:,1:end-1,:) );
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
    ws = fillmask(ws, 0, masks, nnel);
    ud = reshape(A*reshape(us, Np*Mp*Lp,1), Nc,Mc,Lc);
    vd = reshape(A*reshape(vs, Np*Mp*Lp,1), Nc,Mc,Lc);
    wd = reshape(A*reshape(ws, Np*Mp*Lp,1), Nc,Mc,Lc);

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

    %% w back to staggered w points

    w = zeros(Nc+1,Mc,Lc);
%   size(w)
%   size(wd)
    w(2:end-1,:,:) = 0.5*(wd(2:end,:,:)+wd(1:end-1,:,:));
    w(end,:,:) = w(end-1,:,:);



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
      ncwrite(bry_name,'ubar_south',squeeze(ubar(1,:))' ,[1 tout]);
      ncwrite(bry_name,'vbar_south',squeeze(vbar(1,:))' ,[1 tout]);
      ncwrite(bry_name,'zeta_south',squeeze(zetac(1,:))',[1 tout]);
      ncwrite(bry_name,'u_south'   ,squeeze(u(:,1,:))'  ,[1 1 tout]);
      ncwrite(bry_name,'v_south'   ,squeeze(v(:,1,:))'  ,[1 1 tout]);
      ncwrite(bry_name,'w_south'   ,squeeze(w(:,1,:))'  ,[1 1 tout]);
    end
    if (bnd == 2)
      ncwrite(bry_name,'ubar_east',squeeze(ubar(:,Lu))  ,[1 tout]);
      ncwrite(bry_name,'vbar_east',squeeze(vbar(:,Lv))  ,[1 tout]);
      ncwrite(bry_name,'zeta_east',squeeze(zetac(:,Lc)) ,[1 tout]);
      ncwrite(bry_name,'u_east'   ,squeeze(u(:,:,Lu))'  ,[1 1 tout]);
      ncwrite(bry_name,'v_east'   ,squeeze(v(:,:,Lv))'  ,[1 1 tout]);
      ncwrite(bry_name,'w_east'   ,squeeze(w(:,:,end))'  ,[1 1 tout]);
    end
    if (bnd == 3)
      ncwrite(bry_name,'ubar_north',squeeze(ubar(Mu,:))' ,[1 tout]);
      ncwrite(bry_name,'vbar_north',squeeze(vbar(Mv,:))' ,[1 tout]);
      ncwrite(bry_name,'zeta_north',squeeze(zetac(Mc,:))',[1 tout]);
      ncwrite(bry_name,'u_north'   ,squeeze(u(:,Mu,:))'  ,[1 1 tout]);
      ncwrite(bry_name,'v_north'   ,squeeze(v(:,Mv,:))'  ,[1 1 tout]);
      ncwrite(bry_name,'w_north'   ,squeeze(w(:,end,:))'  ,[1 1 tout]);
    end
    if (bnd == 4)
      ncwrite(bry_name,'ubar_west',squeeze(ubar(:,1))  ,[1 tout]);
      ncwrite(bry_name,'vbar_west',squeeze(vbar(:,1))  ,[1 tout]);
      ncwrite(bry_name,'zeta_west',squeeze(zetac(:,1)) ,[1 tout]);
      ncwrite(bry_name,'u_west'   ,squeeze(u(:,:,1))'  ,[1 1 tout]);
      ncwrite(bry_name,'v_west'   ,squeeze(v(:,:,1))'  ,[1 1 tout]);
      ncwrite(bry_name,'w_west'   ,squeeze(w(:,:,1))'  ,[1 1 tout]);
    end

  end    % End loop on par_tind

 end    % End loop bnd

 fclose('all')

  return























