%
%  This script takes hraw from the child grid file and applies log-based smoothing
%  to satisfy the max r-factor criterium (rmax) and sets a minimum depth (hmin)
%
%  Afterwards it matches child topography with that of the parent on each of its 
%  boundaries taking into account land masked areas.
%
%-------------------------------------------------------------------


% ROMS parent and child grid directories
  pdir = '/home/nmolem/USW4/';
  cdir = '/home/nmolem/SMCAL/';
  pdir = '/batavia/nmolem/GOMRE/';
  cdir = '/batavia/nmolem/GOMRE/';
% ROMS parent and child grid files
  pgrid = 'USW4_grd.nc';
  cgrid = 'SMCAL_grd.nc';
  pgrid = 'gom5km-subgrid.nc';
  cgrid = 'gom500_grd.nc';

% When matching boundary topo: Only match to parent topography on open boundaries
  obcflag = [1 1 1 1];      % open boundaries flag (1=open , [S E N W])

  pgrid = [pdir pgrid] 
  cgrid = [cdir cgrid] 

   rmax = 0.2;
   hmin = 20;

   gridfile = cgrid;
   lsmooth

   if exist('pgrid')
    mod_cgrid
   end


 

   
   
