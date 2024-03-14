


% ROMS parent and child grid directories
  pdir = '/zulu/nmolem/NEPAC/';
  cdir = '/zulu/nmolem/SMB/';
  pgrid = 'nepac_grd.nc';
  cgrid = 'smb_grd.nc';

  pgrid = [pdir pgrid] 
  cgrid = [cdir cgrid] 


   rmax = 0.2;
   hmin = 5;
   offset = 0.0;

   gridfile = cgrid;
   lsmooth

   if exist('pgrid')
%  Match boundary topo: Only match to parent topography on open boundaries
    obcflag              = [1 1 0 1];      % open boundaries flag (1=open , [S E N W])
    mod_cgrid
   end


 

   
   
