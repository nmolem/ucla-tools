function  [h_bnd,bdimx] = get_hbnd(vname,h);
  % Get the proper bathymetry for variable and boundary

    if contains(vname,'south')
      bdimx = 'xi_rho';
      if contains(vname,'u_')||contains(vname,'ubar_')||contains(vname,'up_')
        h_bnd = 0.5*(h(2:end,1) + h(1:end-1,1));
	bdimx = 'xi_u';
      elseif contains(vname,'v_')
        h_bnd = 0.5*(h(:,1) + h(:,2));
      else
        h_bnd = h(:,1);
      end
    elseif contains(vname,'north')
      bdimx = 'xi_rho';
      if contains(vname,'u_')||contains(vname,'ubar_')||contains(vname,'up_')
        h_bnd = 0.5*(h(2:end,end) + h(1:end-1,end));
	bdimx = 'xi_u';
      elseif contains(vname,'v_')
        h_bnd = 0.5*(h(:,end) + h(:,end-1));
      else
        h_bnd = h(:,end);
      end
    elseif contains(vname,'west')
      bdimx = 'eta_rho';
      if contains(vname,'u_')
        h_bnd = 0.5*(h(1,:) + h(2,:));
      elseif contains(vname,'v_')||contains(vname,'vbar_')||contains(vname,'vp_')
        h_bnd = 0.5*(h(1,2:end) + h(1,1:end-1));
        bdimx = 'eta_v';
      else
        h_bnd = h(1,:);
      end
    elseif contains(vname,'east')
      bdimx = 'eta_rho';
      if contains(vname,'u_')
        h_bnd = 0.5*(h(end,:) + h(end-1,:));
      elseif contains(vname,'v_')||contains(vname,'vbar_')||contains(vname,'vp_')
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
