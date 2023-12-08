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
    par_dir    = '/batavia/NWATL/';
    par_grd    = 'nwatl_grd.nc';
    par_file   = 'nwatl';
    par_thetas = 6.5;
    par_thetab = 1.5;
    par_hc     = 300.0;
    par_N      = 42;
    par_scoord = 'new2008';   % parent type scoord

    % Child...
    chd_dir    = '/batavia/GULF3/';
    chd_grd    = 'gulf25_grd.nc';
    chd_thetas = 6.5;
    chd_thetab = 1.5;
    chd_hc     = 250.0;
    chd_N      = 42;
    chd_scoord = 'new2008';    % child 'new' or 'old' type scoord

%---------------------------------------------------------------------------------------
% 2. BOUNDARY FILE
%---------------------------------------------------------------------------------------
    obcflag              = [1 1 0 0];      % open boundaries flag (1=open , [S E N W])
    bry_cycle            =  0;             % 0 means no cycle
    bry_filename         = 'gulf25_bry.nc';  % bry filename
    bry_type             = 'his';          % 'avg', 'his' or 'rst'
    total_num_files      =  24;            % number of files to read
    first_file           =  360;            % first avg/his file, eg roms_avg.0282.nc gives 282
    first_record         =  1;             % desired record no. from first avg/his file
    last_record          =  6;             % desired record no. from last avg/his file
    num_records_per_file =  6;             % number of records per parent output file

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
    par_data_path = [par_dir par_file '_' bry_type '.'];

    % Create the bry file
    disp(['Creating boundary file: ' bry_filename]);
    r2r_create_bry(bry_filename,childgrid,obcflag,chdscd,bry_cycle);

    % Get parent subgrid bounds
    disp(' ')
    disp('Get parent subgrids for each open boundary')
    limits = r2r_subgrid(parentgrid,childgrid,obcflag) 
    return

    chd_tind = [0];
    for parfile = 1:total_num_files; % loop through all parent data files
       fld  = (parfile-1) * num_records_per_file + first_file;
       sfld = sprintf('%04d',fld);
       par_data = [par_data_path sfld '.nc'];

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
       chd_tind = par_tind - par_tind(1) + 1 + chd_tind(end);

       display('==============================================================')
       display(' ')
       display(['Processing: ' par_data])
       disp(' ')
       display(['--- using parent record(s)  ', num2str(par_tind)])
       display(['--- setting child record(s) ', num2str(chd_tind)]);

       r2r_hv(parentgrid, par_data, par_tind,           ...
              childgrid, bry_filename, chd_tind,     ...
              chdscd,parscd,obcflag,limits)

    end   %  loop through parent data files
    disp(' ')
    disp('=============== Boundary file done ===============')
    disp(' ')
