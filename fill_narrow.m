gname = 'pachug_grd.nc';
mask = ncread(gname,'mask_rho');
mask(mask<0) = 0;

mask_old = mask;


for it = 1:10
fill = mask;
fill(2:end  ,:) = fill(2:end  ,:) + mask(1:end-1,:);
fill(1:end-1,:) = fill(1:end-1,:) + mask(2:end  ,:);
fill(mask<1) = 0;

nf = sum(fill(:)==1);
disp(['filling: ',num2str(nf),' points in 1p NS passages'])
mask(fill==1) = 0;
if nf==0 ;break;end

fill = mask;
fill(:,2:end  ) = fill(:,2:end  ) + mask(:,1:end-1);
fill(:,1:end-1) = fill(:,1:end-1) + mask(:,2:end);
fill(mask<1) = 0;

nf = sum(fill(:)==1);
disp(['filling: ',num2str(nf),' points in 1p ES passages'])
mask(fill==1) = 0;
if nf==0 ;break;end
end

disp('Filling holes')
[ny,nx] = size(mask);

reg = bwlabel(mask,4);

lint =  0; %% size of largest region
lreg = 0; %% number of largest region
nreg = max(max(reg)); %% number of regions
for i = 1:nreg
  int = sum(sum(reg==i));
  if int>lint
    lreg = i;
    lint = int;
  end
end


for ireg = 1:nreg
  if ireg~=lreg
    %% check before setting to zero
    int = sum(sum(reg==ireg));
    if int>nx*ny/10
      disp(['region: ' num2str(ireg) 'is large.'])
    else
     mask(reg==ireg) = 0;
    end
  end
end

ncwrite(gname,'mask_rho',mask);


umask = mask(1:end-1,:).*mask(2:end,:);
umskr = u2rho(umask);
ufill = umask;
ufill(:,2:end  ) = ufill(:,2:end  ) + umask(:,1:end-1);
ufill(:,1:end-1) = ufill(:,1:end-1) + umask(:,2:end);
ufill(umask<1) = 0;



return

