%   Based on the ideas of Sasha Shchepetkin
%   (c) Jeroen Molemaker UCLA, 2008

      
   if 0
     gridfile = '/batavia/SOLOMON/solom_grd.nc';
     rmax = 0.2;
     hmin = 20;
   end

   nc = netcdf(gridfile,'w');
   h = nc{'hraw'}(:);

      if (rmax>0.D0)
        rmax_log=log((1.D0+rmax*0.9)/(1.D0-rmax*0.9));
      else
        rmax_log=0.D0;
      end

      h(h<hmin) = hmin;

      hl = log(h/hmin);
      hl(hl<0) = 0;

      cf1 = double(1.)/6;
      cf2 = 0.25;

      rt = [];
      disp('iter   max(r)')
      iter_max = 200;
      for iter=1:iter_max

        cff = hl(2:end,2:end-1) - hl(1:end-1,2:end-1);
        cr  = abs(cff);
        Op1 = 1.0*cff.*(1-rmax_log./cr);
        Op1(cr < rmax_log) = 0;

        cff = hl(2:end-1,2:end) - hl(2:end-1,1:end-1);
        cr  = abs(cff);
        Op2 = 1.0*cff.*(1-rmax_log./cr);
        Op2(cr < rmax_log) = 0;

        cff = hl(2:end,2:end) - hl(1:end-1,1:end-1);
        cr  = abs(cff);
        Op3 = 1.0*cff.*(1-rmax_log./cr);
        Op3(cr < rmax_log) = 0;

        cff = hl(1:end-1,2:end) - hl(2:end,1:end-1);
        cr  = abs(cff);
        Op4 = 1.0*cff.*(1-rmax_log./cr);
        Op4(cr < rmax_log) = 0;

        hl(2:end-1,2:end-1) = hl(2:end-1,2:end-1) + cf1*( ...
                              Op1(2:end,:)    -Op1(1:end-1,:)      +Op2(:,2:end)      -Op2(:,1:end-1) +   ...
                         cf2*(Op3(2:end,2:end)-Op3(1:end-1,1:end-1)+Op4(1:end-1,2:end)-Op4(2:end,1:end-1) ) );

%% No gradient at the domain boundaries, this is required for consistency with the ROMS open boundary conditions
        hl(1,:)  = hl(2,:);
        hl(end,:)= hl(end-1,:);

        hl(:,1)  = hl(:,2);
        hl(:,end)= hl(:,end-1);

        h=hmin*exp(hl);

%%  Just add constraints to h here: However, if those constraints imply rt>rmax
%%  you will not exit this loop until the max iterations.

%       h(614:616,268:269) = hraw(614:616,268:269);
%       h(615,268) = hraw(615,268);
%       hl(615,268) = log(h(615,268)/hmin);
%       h(613,269) = h(613,269);

        hr = h;
        rfact;
         rt = [rt 2*max(max(r))];
        disp([' ' num2str(iter) '   '  num2str(rt(end))])
        if rt(end)<rmax
          break
        end

      end 

      disp(['writing lsmoothed h to: ' gridfile])
      nc{'h'}(:) = h;
      close(nc)
