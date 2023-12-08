
  % It's possible that, due to small mismatches between the 
  % masks of the parent and the child, some values in the bry
  % file are nan (fully inside the parent mask) while not masked
  % in the child. We need therefore a fill of the nans in the
  % boundary file

  fname = 'nepac_bry_Y2013M01.nc';

  inf = ncinfo(fname);

  nvar = length(inf.Variables)

  for i = 2:nvar
    vname = inf.Variables(i).Name
    f = ncread(fname,vname);
    fi= f;
    x = [1:size(f,1)];
    if ndims(f)>2 % (x,z,t)
     mask = isnan(f(:,end,100))|abs(f(:,end,100))>1e10;
    else % (x,t)
     mask = isnan(f(:,100))|abs(f(:,100))>1e10;
    end
    xm = x; xm(mask)= [];
    if ndims(f)>2 % (x,z,t)
      for k = 1:size(f,2)
       for it = 1:size(f,3)
        fm = f(:,k,it); fm(mask) = [];
        fi(:,k,it) = interp1(xm,fm,x,'nearest','extrap');
       end
      end
    else % (x,t)
      for it = 1:size(f,2)
       fm = f(:,it); fm(mask) = [];
       fi(:,it) = interp1(xm,fm,x,'nearest','extrap');
      end
    end
    ncwrite(fname,vname,fi);
  end
