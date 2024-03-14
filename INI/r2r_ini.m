%---------------------------------------------------------------------------------------
%
%  roms2roms_ini
%
%  Generate initial file from ROMS parent data. 
%
%---------------------------------------------------------------------------------------
% USER-DEFINED VARIABLES & OPTIONS START HERE
%---------------------------------------------------------------------------------------
%

    % Parent...
      par_dir    = '/zulu/nmolem/PACHUG/';
      par_grd    = 'pachug_grd.nc';
      par_file   = 'pachug_rst.20161231235500.nc';
      par_tind   = 1;            % frame number in parent file
      par_thetas = 6.0;
      par_thetab = 6.0;
      par_hc     = 250.0;
      par_N      = 100;
      parscoord  = 'new2008';    % parent 'new' or 'old' type scoord

    % Child...
      chd_dir    = '/zulu/nmolem/NEPAC/';
      chd_grd    = 'nepac_grd.nc';
      chd_file   = 'nepac_ini.20161231.nc';        % name of new ini file
      chd_thetas = 6.0;
      chd_thetab = 6.0;
      chd_hc     = 250.0;
      chd_N      = 100;
      chdscoord  = 'new2008';                 % child 'new' or 'old' type scoord


%---------------------------------------------------------------------------------------
% USER-DEFINED VARIABLES & OPTIONS END HERE
%---------------------------------------------------------------------------------------
%
     % Child and parent s-coord parameters into chdscd and parscd
       chdscd.N       = chd_N;
       chdscd.theta_s = chd_thetas;
       chdscd.theta_b = chd_thetab;
       chdscd.hc      = chd_hc;

       parscd.N       = par_N;
       parscd.theta_s = par_thetas;
       parscd.theta_b = par_thetab;
       parscd.hc      = par_hc;
       parscd.tind    = par_tind;
       parscd.file    = [par_dir par_file];

     % ROMS parent and child grid files
       pargrd = [par_dir par_grd] 
       chdgrd = [chd_dir chd_grd] 

       parini = [par_dir par_file];
       chdini = [chd_dir chd_file];

     if ~exist(chdini)
      disp(['Creating initial file: ' chdini]);
      r2r_create_ini(chdini, chdgrd, chd_N, chdscd)
     end


     par_grd  = pargrd;
     chd_grd  = chdgrd;
     par_data = parini;
     chd_data = chdini;
     r2r_make_ini

     csr = ncreadatt(parini,'/','Cs_r');
     ncwriteatt(chdini,'/','Cs_r',csr);

     csw = ncreadatt(parini,'/','Cs_w');
     ncwriteatt(chdini,'/','Cs_w',csw);

%    r2r_make_ini(pargrd, parini, chdgrd, chdini, chdscd,parscd,parscoord,chdscoord)
