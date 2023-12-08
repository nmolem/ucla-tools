% parent x: [0,266.4 km];   x0 = 133.2 km
% parent y: [0,319.68 km];  y0 =  99.9 km

pname = 'sm_dx370.nc';
cname = 'sm_dx100.nc';

dx = 100;
pxsize = 2.6667e+05;
cxsize = dx*1024;
x0 = (pxsize - dx*1024)/2;

y0 = 2.5e4;

pname = 'sm_dx370.nc';
cname = 'sm_dx100.nc';

x = [-0.5:1:1024+0.5]*dx + x0;
y = [-0.5:1:2048+0.5]*dx + y0;
[x,y] = meshgrid(x,y);

f   = 0*x+ 7e-5;
ang = 0*x;

ncwrite(cname,'lon_rho',x')
ncwrite(cname,'lat_rho',y')
ncwrite(cname,'f',f')

px = ncread(pname,'lon_rho')';
py = ncread(pname,'lat_rho')';
ph = ncread(pname,'h')';

h = interp2(px,py,ph,x,y);
ncwrite(cname,'h',h')
ncwrite(cname,'hraw',h')



