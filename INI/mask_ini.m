

ncg = netcdf('/batavia/GULFS/gulfs_grd.nc');
mask = ncg{'mask_rho'}(:);
close(ncg)
nci = netcdf('/batavia/GULFS/gulfs_ini.nc','w');

umask = mask(:,1:end-1).* mask(:,2:end);
vmask = mask(1:end-1,:).* mask(2:end,:);

nci{'zeta'}(:) = nci{'zeta'}(:).*mask;
nci{'ubar'}(:) = nci{'ubar'}(:).*umask;
nci{'vbar'}(:) = nci{'vbar'}(:).*vmask;
for k = 1:42
	k
 nci{'u'}(1,k,:,:) = squeeze(nci{'u'}(1,k,:,:)).*umask;
 nci{'v'}(1,k,:,:) = squeeze(nci{'v'}(1,k,:,:)).*vmask;
 nci{'temp'}(1,k,:,:) = squeeze(nci{'temp'}(1,k,:,:)).*mask;
 nci{'salt'}(1,k,:,:) = squeeze(nci{'salt'}(1,k,:,:)).*mask;
end
close(nci)

