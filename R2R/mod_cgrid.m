%--------------------------------------------------------------
%
%  Modify child grid topography such that it matches the interpolated
%  parent topography at the boundaries.
%
%  This script is for use with make_roms2roms.
%
%   (c) 2007 Jeroen Molemaker, UCLA
%   (c) 2014 Updated for intrinsic matlab netcdf functionality Jeroen Molemaker, UCLA
%--------------------------------------------------------------
%
if 0
  clear all
% ROMS parent and child grid directories
  pdir = '/home/nmolem/USW4/';
  cdir = '/home/nmolem/SMCAL/';
% ROMS parent and child grid files
  pgrid = 'usw4_grd.nc';
  cgrid = 'smcal_grd.nc';

% Only match to parent topography on open boundaries
  obcflag              = [1 1 1 1];      % open boundaries flag (1=open , [S E N W])

  pgrid = [pdir pgrid];
  cgrid = [cdir cgrid];
 end
% End user-defined----------------------------------------------
%

% Get minimal parent subgrid bounds
    nc = netcdf(cgrid,'nowrite');
    lonc = nc{'lon_rho'}(:);
    latc = nc{'lat_rho'}(:);
    close(nc);
%   lonc = ncread(cgrid,'lon_rho')';
%   latc = ncread(cgrid,'lat_rho')';

    nc    = netcdf(pgrid,'nowrite');
    lonp  = nc{'lon_rho'}(:);
    latp  = nc{'lat_rho'}(:);
    close(nc);
%   lonp = ncread(pgrid,'lon_rho')';
%   latp = ncread(pgrid,'lat_rho')';

    lon0 = min(min(lonc))-0.05;
    lon1 = max(max(lonc))+0.05;
    lat0 = min(min(latc))-0.05;
    lat1 = max(max(latc))+0.05;

    g = lonp>=lon0&lonp<=lon1 & latp>=lat0&latp<=lat1;
    jmin = min(find(any(g')));
    jmax = max(find(any(g')));
    imin = min(find(any(g)));
    imax = max(find(any(g)));
    clear g

   
    plot(lonp(jmin:jmax,imin:imax),latp(jmin:jmax,imin:imax),'.k')
    hold on;plot(lonc,latc,'.r');hold off
    drawnow

% Get topography data from childgrid
  nc = netcdf(cgrid,'nowrite');
  hc   = nc{'h'}(:);
  lonc = nc{'lon_rho'}(:);
  latc = nc{'lat_rho'}(:);
  mask = nc{'mask_rho'}(:);
  close(nc);
% hc   = ncread(cgrid,'h')';
% mask = ncread(cgrid,'mask_rho')';
  [Mc,Lc]=size(hc);

% Get parent grid and squeeze minimal subgrid
  nc    = netcdf(pgrid,'nowrite');
  hp    = squeeze(nc{'h'}(jmin:jmax,imin:imax));
  lonp  = squeeze(nc{'lon_rho'}(jmin:jmax,imin:imax));
  latp  = squeeze(nc{'lat_rho'}(jmin:jmax,imin:imax));
  maskp = squeeze(nc{'mask_rho'}(jmin:jmax,imin:imax));
  lonp(lonp<0) = lonp(lonp<0);
  close(nc);
% hp   = ncread(pgrid,'h')';
% maskp= ncread(pgrid,'mask_rho')';

% hp    = hp(jmin:jmax,imin:imax);
% lonp  = lonp(jmin:jmax,imin:imax);
% latp  = latp(jmin:jmax,imin:imax);
% maskp = maskp(jmin:jmax,imin:imax);
 [Mp,Lp]=size(hp);

% Get interpolation coefficient to go to (lonc,latc).
  [elem,coef] = get_tri_coef(lonp,latp,lonc,latc,maskp);
%% parent grid topo at child locations
   hpi = sum(coef.*hp(elem),3);


  if 0
  dist = zeros(Mc,Lc,4);
  for i = 1:Mc   %% north south
   for j = 1:Lc    %% east west
     dist(i,j,1) =      i/Mc + (1-obcflag(1))*1e6; % South
     dist(i,j,2) = (Lc-j)/Lc + (1-obcflag(2))*1e6; % East
     dist(i,j,3) = (Mc-i)/Mc + (1-obcflag(3))*1e6; % North
     dist(i,j,4) =      j/Lc + (1-obcflag(4))*1e6; % West
   end
  end
  dist = min(dist,[],3);

  alpha = 0.5*tanh(100*(dist-0.03))+0.5; %% Feel free to play with this function.
  alpha = 0.5*tanh( 50*(dist-0.06))+0.5; %% Feel free to play with this function.
  else
   distance
  end

  hcn = alpha.*hc + (1-alpha).*hpi;

if 1
  nc = netcdf(cgrid,'write');
  nc{'h'}(:) = hcn;
  nc{'h'}.notes1 = 'Topo has been modified to match the parent grid topo at the boundaries';
  close(nc);
 end

%% Visualize the modification
  sc0 = min(min(hcn));
  sc1 = max(max(hcn));
  subplot(2,2,1)
  pcolor(lonc,latc,hpi);caxis([sc0 sc1]);colorbar;shading flat
  title('Interpolated Parent Topo')
  subplot(2,2,2)
  pcolor(lonc,latc,hcn);caxis([sc0 sc1]);colorbar;shading flat
  title('Boundary Smoothed Child Topo')
  subplot(2,2,3)
  pcolor(lonc,latc,hcn-hpi);colorbar;shading flat
  title('Difference between Parent and child Topo');
  subplot(2,2,4)
  pcolor(lonc,latc,alpha);colorbar;shading flat
  title('Parent/Child transition function');
