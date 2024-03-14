   gridfile = '/zulu/nmolem/GWAKE/gwake_grd.nc';
%  Make 
%  clear jet
%  colormap(jet(256))

   UTM=0;

   d2r = pi/180;
   r2d = 180/pi;
   r_earth = 6371315.;

%  1. Extract minimal set from topo file.  (plus 'width' size borders)
%  2. Establish bounding rectangle.
%  3. Based on weighting scale, determine maximum stencil.

%  TODO: add wrapper border whenever the grid stradles the dateline.
   % input parameters
   formask = 0;
   if formask
     wd = 1.;  % for mask
     climat = 'taidp200';
   else
     wd= 32;  % for hraw
     climat = 'taidp200';
   end
   wd = wd/2;

   if strcmp(lower(climat),'sfbay_10m')||strcmp(lower(climat),'sfbay_cen_1m')
    UTM=1;
   end

   global_grid = 1;
   if global_grid
    xg_glob = ncread(gridfile,'lon_rho');
    yg_glob = ncread(gridfile,'lat_rho');
    pn_glob = ncread(gridfile,'pn');
    h_glob = xg_glob + nan;
   end
   if UTM
    xg_glob = xg_glob-360;
    [xg_glob,yg_glob] = ll2utm(xg_glob,yg_glob);
    xg_glob = xg_glob*1e3;
    yg_glob = yg_glob*1e3;
   end

   [nxg,nyg] = size(xg_glob);

   display(['Smoothing radius: ' num2str(wd) ' gridpoints'])

   if strcmp(lower(climat),'taidp200') & ~formask
    ncwriteatt(gridfile,'hraw','source',['Raw bathymetry from Taidp200m (smoothing radius ' num2str(wd) ')']);
   end
   if strcmp(lower(climat),'sfbay_cen_1m')  & ~formask
    ncwriteatt(gridfile,'hraw','source',['Raw bathymetry from 1m DEM for Cental SF bay (www.sciencebase.gov) (Smoothing radius ' num2str(wd) ' )']);
   end


   nchx =  1;
   nchy =  1;
   chx = round([0:nchx]*nxg/nchx);
   chy = round([0:nchy]*nyg/nchy);
   ci = 1;
   cj = 1;

  xl0 = 10000000;
  xl1 =-10000000;
  yl0 = 10000000;
  yl1 =-10000000;
% for  ci = 3:3
%  for  cj = 2:2
  for  ci = 1:nchx
   for  cj = 1:nchy
     tic
     disp(['chunk ' num2str(ci) num2str(cj)])
     gj0 = chy(cj)+1;
     gj1 = chy(cj+1);
     gi0 = chx(ci)+1;
     gi1 = chx(ci+1);

%    ncg = netcdf(gridfile,'w');
     xg = xg_glob(gi0:gi1,gj0:gj1);
     yg = yg_glob(gi0:gi1,gj0:gj1);
     pn = pn_glob(gi0:gi1,gj0:gj1);
     [nxg,nyg] = size(xg);

     dg = 1./min(min(pn));
     width = wd*dg;  %% largest smoothing width for this patch

     %% figure out how large the patch of topo data needs to be
     xg_min = min(min(xg));
     xg_max = max(max(xg));
     yg_min = min(min(yg));
     yg_max = max(max(yg));

     coslat = cos(d2r*max(max(abs(yg))));
     dellon = r2d*width/(coslat*r_earth);
     dellat = r2d*width/(r_earth);

     if UTM
      dellon = width;
      dellat = width;
     end 

     xg_min = xg_min-1.2*dellon*2;
     xg_max = xg_max+1.2*dellon*2;
     yg_min = yg_min-1.2*dellat*2;
     yg_max = yg_max+1.2*dellat*2;

     %% Just for plotting
%    if UTM
%     xl0 = xg_min;
%     xl1 = xg_max;
%     yl0 = yg_min;
%     yl1 = yg_max;
%    else
      xl0 = min(xl0,xg_min);
      xl1 = max(xl1,xg_max);
      yl0 = min(yl0,yg_min);
      yl1 = max(yl1,yg_max);
%    end

     switch lower(climat);
      case{'sfbay_cen_1m'}
        disp('Using SFbay_Cen 1m dem')
        datafile = '/zulu/nmolem/DATASETS/TOPO/SFbay_Cen_1m_dem.nc'; %% from 
        topo_lon = ncread(datafile,'x'); 
        topo_lat = ncread(datafile,'y'); 
      case{'sfbay_10m'}
        disp('Using SFbay 10m dem')
        datafile = '/zulu/nmolem/DATASETS/TOPO/SFbay_10m_dem.nc'; %% from 
        topo_lon = ncread(datafile,'x');
        topo_lat = ncread(datafile,'y');
      case{'taidp200'}
        disp('Using Taidp200')
        datafile = '/zulu/nmolem/DATASETS/TOPO/Taidp200.nc'; %% from 
        topo_lon = ncread(datafile,'longitude');
        topo_lat = ncread(datafile,'latitude');
      case{'taidp500'}
        disp('Using Taidp500')
        datafile = '/zulu/nmolem/DATASETS/TOPO/Taidp500.nc'; %% from 
        topo_lon = ncread(datafile,'longitude');
        topo_lat = ncread(datafile,'latitude');
      case{'crm_west'}
        disp('Using CRM_west')
        datafile = '/batavia/nmolem/OBSERV/TOPO/crm_west.nc'; %% from 
        topo_lon = ncread(datafile,'x');
        topo_lat = ncread(datafile,'y');
        topo_lon = topo_lon + 360;
      case{'usgs_chesa'}
        disp('Using USGS data set')
        datafile = '/avatar/nmolem/batavia/nmolem/OBSERV/TOPO/usgs_chesapeake.nc'; %% from -77 to -74 degrees
        topo_lon = ncread(datafile,'lon');
        topo_lat = ncread(datafile,'lat');
      case{'etopo'}
        disp('Using Etopo1 data set')
%       nct = netcdf('/batavia/nmolem/OBSERV/TOPO/ETOPO2v2c_f4.nc'); %% from -180 180 degrees
%       nct = netcdf('/batavia/nmolem/OBSERV/TOPO/etopo2_0_360.nc'); %% from 0 360 degrees
%       topo_lon = nct{'x'}(:);
%       topo_lat = nct{'y'}(:);

        datafile = '/avatar/nmolem/batavia/nmolem/OBSERV/TOPO/ETOPO1_180.nc'; %% from -180 180 degrees
        topo_lon = ncread(datafile,'lon');
        topo_lat = ncread(datafile,'lat');
      case{'gebco'}
        disp('Using Gebco data set')
        nct = netcdf('/batavia/nmolem/OBSERV/TOPO/topo_gebco_GRANCAN.nc');
        topo_lon = nct{'lon'}(:);
        topo_lat = nct{'lat'}(:);
      case{'srtm30'}
        disp('Using SCRIPS data set')
%       datafile = '/avatar/nmolem/batavia/nmolem/OBSERV/TOPO/SRTM30_180.nc';
        datafile = '/avatar/nmolem/batavia/nmolem/OBSERV/TOPO/SRTM30_0_360_new.nc';
        topo_lon = ncread(datafile,'longitude');%   - 360;
        topo_lat = ncread(datafile,'latitude');
      case{'palau'}
        datafile = '/avatar/nmolem/batavia/nmolem/OBSERV/TOPO/palau_bathy.nc';
        topo_lon = ncread(datafile,'lon');
        topo_lat = ncread(datafile,'lat');
	topo_lon = topo_lon(:,1);
	topo_lat = topo_lat(1,:);
      case{'srtm15'}
        disp('Using the SRTM15 data set')
        datafile = '/avatar/nmolem/batavia/nmolem/OBSERV/TOPO/SRTM15_V2.5.nc';
        topo_lon = ncread(datafile,'lon');
        topo_lat = ncread(datafile,'lat');
        topo_lon = topo_lon + 360;
      case{'socal_30m'}
        disp('Using the Socal 30m data set')
        datafile = '/batavia/nmolem/OBSERV/TOPO/socal_30m.nc';
        topo_lon = ncread(datafile,'lon');
        topo_lat = ncread(datafile,'lat');
        topo_lon = topo_lon + 360;
     end

     dx_topo = r_earth*(topo_lon(2)-topo_lon(1))*d2r;
     dy_topo = r_earth*(topo_lat(2)-topo_lat(1))*d2r;
     if UTM
      dx_topo = topo_lon(2)-topo_lon(1);
      dy_topo = topo_lat(2)-topo_lat(1);
     end

%     [min(topo_lon) xg_min]
     if max(topo_lat)<yg_min|min(topo_lat)>yg_max
       disp('chunk not in range')
       continue
     end
     if max(topo_lon)<xg_min|min(topo_lon)>xg_max
       disp('chunk not in range')
       continue
     end

     if ~strcmp(lower(climat),'sands')
      ibt = find( topo_lon > xg_min,1,'first')-1;
      iet = find( topo_lon < xg_max,1,'last')+1;
      jbt = find( topo_lat > yg_min,1,'first')-1;
      jet = find( topo_lat < yg_max,1,'last')+1;
      ibt = max(ibt,1); iet = min(iet,length(topo_lon));
      jbt = max(jbt,1); jet = min(jet,length(topo_lat));
%     if yg_min<min(topo_lat)+0.04
%      jbt = 1;
%      yg(yg<min(topo_lat)+0.04) = min(topo_lat)+ 0.04;
%     end
%     if xg_min<min(topo_lat)+0.04
%      ibt = 1;
%      xg(xg<min(topo_lon)+0.04) = min(topo_lon)+ 0.04;
%     end
      topo_lon = topo_lon(ibt:iet);
      topo_lat = topo_lat(jbt:jet);
     end

  switch lower(climat);
    case{'sfbay_cen_1m'}
      hin = -ncread(datafile,'z',[ibt jbt],[iet-ibt+1 jet-jbt+1]);
    case{'sfbay_10m'}
      hin = -ncread(datafile,'z',[ibt jbt],[iet-ibt+1 jet-jbt+1]);
    case{'taidp200'}
      hin = -ncread(datafile,'elevation',[ibt jbt],[iet-ibt+1 jet-jbt+1]);
    case{'taidp500'}
      hin = -ncread(datafile,'elevation',[ibt jbt],[iet-ibt+1 jet-jbt+1]);
    case{'crm_west'}
%     hin = nct{'topo'}(jbt:jet,ibt:iet);
      hin = -ncread(datafile,'z',[ibt jbt],[iet-ibt+1 jet-jbt+1]);
    case{'usgs_chesa'}
%     hin = nct{'topo'}(jbt:jet,ibt:iet);
      hin = ncread(datafile,'topo',[ibt jbt],[iet-ibt+1 jet-jbt+1]);
    case{'etopo'}
      hin = -ncread(datafile,'z',[ibt jbt],[iet-ibt+1 jet-jbt+1]);
    case{'gebco'}
      hin = -nct{'topo'}(jbt:jet,ibt:iet);
    case{'srtm30'}
%     hin = -nct{'elevation'}(jbt:jet,ibt:iet);
      hin = -ncread(datafile,'elevation',[ibt jbt],[iet-ibt+1 jet-jbt+1]);
%     hin(hin>7500) = 7500;disp('limiting srtm30 topo at 7500m');
    case{'sands'}
      [hin,topo_lat,topo_lon] = extract_1m([yg_min-0.02 yg_max+0.02 xg_min-0.02 xg_max+0.02],1);
       topo_lat = topo_lat(1:end-1);
       hin = -hin(2:end,:);
    case{'palau'}
      hin = ncread(datafile,'fdept',[ibt jbt],[iet-ibt+1 jet-jbt+1]);
    case{'srtm15'}
      hin = -ncread(datafile,'z',[ibt jbt],[iet-ibt+1 jet-jbt+1]);
    case{'socal_30m'}
      hin = -ncread(datafile,'z',[ibt jbt],[iet-ibt+1 jet-jbt+1]);
   end
   [nxt,nyt] = size(hin);

   if ~formask
%   hin(hin<5) = 5;
%   hin(hin<5) = 5;
   end

   if sum(~isnan(hin(:)))<10
     display('no data for this chunk')
     continue
   end

   if 0  %% testing
     [xt2,yt2] = meshgrid(xt,yt);
     hin = xt2;
   end

   ill = xg;
   jll = yg;
   lon_min = min(topo_lon);
   lat_min = min(topo_lat);
   for i = 1: nxg  %% Locate a index in the topo file that is close to the grid point
    for j = 1: nyg
      if lat_min>yg(i,j) 
        jll(i,j) = 1;
      else
        jll(i,j) = find( topo_lat < yg(i,j),1,'last');
      end
      if lon_min>xg(i,j)
        ill(i,j) = 1;
      else
        ill(i,j) = find( topo_lon < xg(i,j),1,'last');
      end
    end 
   end 
   sl = size(topo_lat);


   %% minimal distance of dx is:   
   %%  xt and yt should be monotonous and increasing
    lcoslat = cos(d2r*max(abs(topo_lat)));

    if UTM
     isz = 2*ceil(width/dx_topo)
     jsz = 2*ceil(width/dy_topo)
    else
     isz = 2*ceil(width/dx_topo/lcoslat) 
     jsz = 2*ceil(width/dy_topo)
    end
%   siz = max(isz,jsz);
%   isz = siz;jsz=siz;

 if isz<3 & jsz<3
   display('Doing linear interpolation');
   [xt,yt] = meshgrid(topo_lon,topo_lat);
%  xt = xt';
%  yt = yt';
%  hraw = interp2(xt,yt,hin,xg,yg,'makima');
   hraw = interp2(xt,yt,hin',xg,yg,'linear');
 else
   display('Doing Weighted average');

   shfti = isz/2 -1;
   shftj = jsz/2 -1;
   weight = zeros(isz,jsz);
   hraw   = zeros(nxg,nyg) + nan;
   for i = 1:nxg
    for j = 1:nyg
      ib = ill(i,j)      -shfti;
      ie = ill(i,j)+isz-1-shfti;

      jb = jll(i,j)      -shftj;
      je = jll(i,j)+jsz-1-shftj;

      ib = max(ib,1);ie = min(ie,nxt);
      jb = max(jb,1);je = min(je,nyt);
       
      xtl = topo_lon(ib:ie);
      ytl = topo_lat(jb:je);
      [xtl,ytl] = meshgrid(xtl,ytl);
      xtl = xtl'; ytl = ytl';
      if UTM
       r = sqrt( (xtl-xg(i,j)).^2 + (ytl-yg(i,j)).^2 );
      else
       r = gc_dist(xg(i,j),yg(i,j),xtl,ytl,1);
      end
      lwidth = wd./pn(i,j);  %% local cell size
%     [width lwidth]
      w = (1-(r/lwidth).^2).^2;
      w(r>lwidth) = 0;

      hlocal = hin(ib:ie,jb:je);

%     if 0 % used for crm
%     w(isnan(hlocal)) = 0;
%     hlocal(isnan(hlocal)) = 0;
%     ratio = sum(w(:)==0)/length(w(:));
%     if ratio>0.5
%       continue
%     end
%     end

      w = w/sum(sum(w));
%     return
%     if isnan(w)>0
%     end

      if min([min(w(1,:)) min(w(end,:)) min(w(:,1)) min(w(:,end))]) >0
%       figure(3);contour(w,[0 0],'k');drawnow
%       error 'topo sampling region too small'
        hraw(i,j) = nan;
      else
       hraw(i,j) = sum(w(:).*hlocal(:));
      end


    end 
    if mod(i,10)==0 
%     hraw(hraw<-0.2) = -0.2;
      figure(1);mypcolor(xg,yg,hraw);
      pause(0.01)
    end
   end 
 end 
 %hraw(hraw<-0.2) = -0.2;

% hraw(hraw==0) = nan;

% hraw(hraw<0) = 0;
  mask = 1 + 0.*hraw;
  mask(hraw<=0.11) = 0;

   if formask
     disp('writing mask')
     ncwrite(gridfile,'mask_rho',mask,[gi0 gj0]);
%    ncg{'mask_rho'}(gi0:gi1,gj0:gj1)= mask;
   else
     disp('writing topo')
     ncwrite(gridfile,'hraw',hraw,[gi0 gj0]);
   end

  if global_grid
   h_glob(gi0:gi1,gj0:gj1) =  hraw;
   figure(2);
   mypcolor(xg_glob,yg_glob,h_glob);xlim([xl0 xl1]);ylim([yl0 yl1])
  end
 toc
end  %% chj loop
end  %% chi loop



