function  [h_bnd,bdimx] = get_hbnd(vname,h);
  % Get the proper bathymetry for variable and boundary

    if strfind(vname,'south')
      bdimx = 'xi_rho';
      if strfind(vname,'u_')
        h_bnd = 0.5*(h(2:end,1) + h(1:end-1,1));
	bdimx = 'xi_u';
      elseif strfind(vname,'v_')
        h_bnd = 0.5*(h(:,1) + h(:,2));
      else
        h_bnd = h(:,1);
      end
    elseif strfind(vname,'north')
      bdimx = 'xi_rho';
      if strfind(vname,'u_')
        h_bnd = 0.5*(h(2:end,end) + h(1:end-1,end));
	bdimx = 'xi_u';
      elseif strfind(vname,'v_')
        h_bnd = 0.5*(h(:,end) + h(:,end-1));
      else
        h_bnd = h(:,end);
      end
    elseif strfind(vname,'west')
      bdimx = 'eta_rho';
      if strfind(vname,'u_')
        h_bnd = 0.5*(h(1,:) + h(2,:));
      elseif strfind(vname,'v_')
        h_bnd = 0.5*(h(1,2:end) + h(1,1:end-1));
        bdimx = 'eta_v';
      else
        h_bnd = h(1,:);
      end
    elseif strfind(vname,'east')
      bdimx = 'eta_rho';
      if strfind(vname,'u_')
        h_bnd = 0.5*(h(end,:) + h(end-1,:));
      elseif strfind(vname,'v_')
        h_bnd = 0.5*(h(end,2:end) + h(end,1:end-1));
        bdimx = 'eta_v';
      else
        h_bnd = h(end,:);
      end
    else
        error('Cant determine boundary')
    end

    h_bnd = h_bnd';

    return
