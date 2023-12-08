%---------------------------------------------------------------------------------------
%
%  make_roms2roms
%
%  Generate boundary perimeter file from ROMS parent data.
%  Designed for compatibility with
%  arbitrary horizontal orientation of the grid system and
%  desired vertical s-coordinate ('new' or 'old').
%
%  Note that when run this script it tests for the presence of a .mat file
%  which contains various interpolation coefficients related to your child
%  and parent grids.  If the .mat file is not there it will calculate the coefficients
%
%
%  Jeroen Molemaker and Evan Mason in 2007-2009 at UCLA
%
%---------------------------------------------------------------------------------------
clear all
close all
disp(' ')
%---------------------------------------------------------------------------------------
% USER-DEFINED VARIABLES & OPTIONS START HERE
%---------------------------------------------------------------------------------------
%
%---------------------------------------------------------------------------------------
% 1.  GENERAL
%---------------------------------------------------------------------------------------

    % Parent...
%   par_dir    = '/mnt/modisk/comanche/nmolem/GULF/';
%   par_grd    = 'gulf25_grd.nc';
%   par_dir    = '/mnt/ladisk/ROMS2km/data/';
    par_dir    = '/mnt/ladisk/IASROMS-2km/Phase4c/';
    par_grd    = '2010040100_da.nc';
    par_file   = '2010';
    par_thetas = 7.0;
    par_thetab = 2.0;
    par_hc     = 300.0;
    par_N      = 30;
    par_scoord = 'new2008';   % parent type scoord

    % Child...
    chd_dir    = '/mnt/batavia/NOAA/NEST/';
    chd_grd    = 'gom_l2_grd.nc';
    chd_thetas = 7.0;
    chd_thetab = 2.0;
    chd_hc     = 300.0;
    chd_N      = 40;
    chd_scoord = 'new2008';    % child 'new' or 'old' type scoord

%---------------------------------------------------------------------------------------
% 2. BOUNDARY FILE
%---------------------------------------------------------------------------------------
    obcflag              = [1 1 0 1];      % open boundaries flag (1=open , [S E N W])
    bry_cycle            =  0;             % 0 means no cycle
    bry_filename         = 'gom_l2b_bry.nc'; % bry filename
    bry_type             = 'dhis';         % 'avg', 'his' or 'rst'
    total_num_files      = 366*12;         % number of files to read
    first_file           =   1;            % first avg/his file, eg roms_avg.0282.nc gives 282
    first_record         =   1;            % desired record no. from first avg/his file
    last_record          =   1;            % desired record no. from last avg/his file
    num_records_per_file =   1;            % number of records per parent output file

%---------------------------------------------------------------------------------------
% USER-DEFINED VARIABLES & OPTIONS END HERE
%---------------------------------------------------------------------------------------
%

% Put the various paths/filenames/variables together

    % Child and parent s-coord parameters into chdscd and parscd
    chdscd.theta_s = chd_thetas;
    chdscd.theta_b = chd_thetab;
    chdscd.hc      = chd_hc;
    chdscd.N       = chd_N;
    chdscd.scoord  = chd_scoord;


    parscd.theta_s = par_thetas;
    parscd.theta_b = par_thetab;
    parscd.hc      = par_hc;
    parscd.N       = par_N;
    parscd.scoord  = par_scoord;

    % ROMS parent and child grid files
    parentgrid = [par_dir  par_grd] 
    childgrid  = [chd_dir  chd_grd] 


%---------------------------------------------------------------------------------------
%   BOUNDARY PERIMETER FILE
%---------------------------------------------------------------------------------------

    !rm *.mat

    bry_filename  = [chd_dir bry_filename];
    if 0 %% normal ROMS
    par_data_path = [par_dir par_file '_' bry_type '.'];
    else  % different
    par_data_path = [par_dir par_file];
    end

    % Create the bry file
%   disp(['Creating boundary file: ' bry_filename]);
%   r2r_create_bry(bry_filename,childgrid,obcflag,chdscd,bry_cycle);

    % Get parent subgrid bounds
    disp(' ')
    disp('Get parent subgrids for each open boundary')
    limits = r2r_bry_subgrid(parentgrid, childgrid,obcflag) 

%   chd_tind = [0];
%   for parfile = 1:total_num_files; % loop through all parent data files
%   for parfile = 1:150
    chd_tind = [4250];
    for parfile = 4251:total_num_files
      if 0 %% normal ROMS nomenclature
        fld  = (parfile-1) * num_records_per_file + first_file;
        sfld = sprintf('%04d',fld);
        par_data = [par_data_path sfld '.nc'];
      else  %% special iasroms file names
	 dnum0 = datenum(2010,04,01);
	 dtime = (parfile-1)/24.0 + dnum0;
	 [yr,month,day,hour,mn,sec] = datevec(dtime);

         smn = sprintf('%02d',month);
         sdy = sprintf('%02d',day);
         shr = sprintf('%02d',hour);
         par_data = [par_data_path smn sdy shr '_da.nc'] 
         parscd.par_file  = par_data;
      end

       if (parfile==1)
        rec0 = first_record;
       else
        rec0 = 1.;
       end
       if (parfile==total_num_files)
        rec1 = last_record;
       else
        rec1 = num_records_per_file;
       end
       par_tind = [rec0:rec1];
       chd_tind = par_tind - par_tind(1) + 1 + chd_tind(end) 

       display('==============================================================')
       display(' ')
       display(['Processing: ' par_data])
       disp(' ')
       display(['--- using parent record(s)  ', num2str(par_tind)])
       display(['--- setting child record(s) ', num2str(chd_tind)]);

       r2r_hv(parentgrid, par_data, par_tind,           ...
              childgrid, bry_filename, chd_tind,     ...
              chdscd,parscd,obcflag,limits)

       fclose('all')
    end   %  loop through parent data files
    disp(' ')
    disp('=============== Boundary file done ===============')
    disp(' ')
