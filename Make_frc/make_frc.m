
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Make a ROMS bulk forcing file using hourly data Atmospheric data.
%  Choice of ERA5 (25km), IFS (9km), or HRRR (3km)  (USA coastal regions only)
%
%  2020-2024, Jeroen Molemaker, Pierre Damien, UCLA
%
%
%%%%%%%%%%%%%%%%%%%%% USER-DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
%
%  frc climatology file names:
frc_root = '/zulu/DATASETS/';

frc_source = 'HRRR'; % ERA5/IFS/HRRR

% Set a date range for the forcing files
start_date = datenum(2017,06,02);
end_date   = datenum(2017,08,30);

grdname  = '/zulu/nmolem/SMODE/smode_grd.nc';
root_name= '/zulu/nmolem/SMODE/smode';

% Only needed in case of wind_dropoff=1
disname  = '/zulu/nmolem/XXX/xxx_cdist.mat';

coarse_frc   = 1; % forcing files at half the resolution of the grid
rad_corr     = 0; % Multiplicative correction of swr and lwr to observations

wind_dropoff = 0; % Spatial field to represent coastal wind dropoff 
add_rivers   = 0; % Adds river runoff as additional precipitation (Obsolete)
%
%%%%%%%%%%%%%%%%%%% END USER-DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%

switch frc_source
  case 'ERA5'
    disp('Using ERA5 atmospheric data')
    dsatt = 'ERA5 (25 km nominal res)';

    frc_dir = [frc_root 'ERA5/'];
    maskname = [frc_dir 'ERA5_mask.nc'];
    rcorname = [frc_dir 'ERA5_rad_cor.nc'];
    fill_frc_ecnwf
  case 'IFS'
    disp('Using IFS atmospheric data')
    dsatt = 'IFS (9 km nominal res)';

    frc_dir = [frc_root 'IFS/'];
    rcorname = [frc_dir 'IFS_rad_cor_core.nc'];
    maskname = [frc_dir 'IFS_mask.nc'];
    fill_frc_ecnwf
  case 'HRRR'
    disp('Using HRRR atmospheric data')
    dsatt = 'HRRR (3 km res)';

    frc_dir = [frc_root 'HRRR/'];
    rcorrname = [frc_dir 'none'];
    maskname  = [frc_dir 'hrrr_land.nc'];

    fill_frc_hrrr
  otherwise
    disp('Unknown atmospheric data source')
end
