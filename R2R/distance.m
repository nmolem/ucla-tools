%
%  This is a utility function for mod_cgrid.
%  It computes distance to the nearest non-masked
%  point on the boundary.
%  
%  The distance field is used to compute a parameter [0,1]
%  which merges smoothly between parent and child topography
%
%  (c) 2010, Jeroen Molemaker, UCLA
%---------------------------------------------------------
  dist = mask*0 + 1e6;
  [ny,nx] = size(mask);
  n = max(nx,ny);

  x = [0:nx-1]/n;
  y = [0:ny-1]/n;
  [x,y] = meshgrid(x,y);

  trans = 0.05;
  width = ceil(n*trans);
  bx = [];
  by = [];
  if obcflag(1)==1  %% South boundary
    bx = x(1,:); bx(mask(1,:)<1) = [];
    by = y(1,:); by(mask(1,:)<1) = [];
    nb = length(bx);
    for i = 1: nb
	   [i nb]
     dtmp = (x(1:width,:)-bx(i) ).^2 + (y(1:width,:)-by(i) ).^2;
     dist(1:width,:) = min(dist(1:width,:),dtmp);
    end
  end
  if obcflag(2)==1  %% East boundary
    bx = x(:,end); bx(mask(:,end)<1) = [];
    by = y(:,end); by(mask(:,end)<1) = [];
    nb = length(bx);
    for i = 1: nb
	   [i nb]
     dtmp = (x(:,nx-width:nx)-bx(i) ).^2 + (y(:,nx-width:nx)-by(i) ).^2;
     dist(:,nx-width:nx) = min(dist(:,nx-width:nx),dtmp);
    end
  end
  if obcflag(3)==1  %% North boundary
    bx = x(end,:); bx(mask(end,:)<1) = [];
    by = y(end,:); by(mask(end,:)<1) = [];
    nb = length(bx);
    for i = 1: nb
	   [i nb]
     dtmp = (x(ny-width:ny,:)-bx(i) ).^2 + (y(ny-width:ny,:)-by(i) ).^2;
     dist(ny-width:ny,:) = min(dist(ny-width:ny,:),dtmp);
    end
  end
  if obcflag(4)==1  %% West boundary
    bx = x(:,1); bx(mask(:,1)<1) = [];
    by = y(:,1); by(mask(:,1)<1) = [];
    nb = length(bx);
    for i = 1: nb
     [i nb]
     dtmp = (x(:,1:width)-bx(i) ).^2 + (y(:,1:width)-by(i) ).^2;
     dist(:,1:width) = min(dist(:,1:width),dtmp);
    end
  end
  dist = sqrt(dist); %/trans;
  dist(dist>trans) = trans;
  dist = dist/trans;
  alpha = 0.5 - 0.5*cos(pi*dist);
% alpha = 




