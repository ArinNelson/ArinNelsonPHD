function status = roms_interp_ini(nc_ini_old,nc_grd_old,nc_ini_new,nc_grd_new)

  % Get info from old file
  info = ncinfo(nc_ini_old);

  % New grid dimensions
  xr_new = ncread(nc_grd_new,'lon_rho');
  xi  = size(xr_new,1);
  eta = size(xr_new,2);
  dim_name  = {info.Dimensions.Name};
  dim_str   = {'xi_rho','xi_u','xi_v','eta_rho','eta_u','eta_v'};
  dim_lngth = [xi,xi-1,xi,eta,eta,eta-1];
  for i=1:numel(dim_str)
    info.Dimensions( strcmp(dim_name,dim_str{i})==1 ).Length = dim_lngth(i);
  end
  clear i;
  
  % Do same for variables
  for i=1:numel(info.Variables)
  
    % Dimensions
    for j=1:numel(info.Variables(i).Dimensions)
      jj = find( strcmp(dim_str,info.Variables(i).Dimensions(j).Name)==1 );
      if(~isempty(jj))
        info.Variables(i).Dimensions(j).Length = info.Dimensions( strcmp(dim_str{jj},dim_name)==1 ).Length;
      end
      info.Variables(i).Size(j) = info.Variables(i).Dimensions(j).Length;
    end

  end
  clear i;
  
  % Generate new file
  if(exist(nc_ini_new,'file')==2); delete(nc_ini_new); end
  ncwriteschema(nc_ini_new,info);
  
  % Transfer non-grid variables
  var_name = {'spherical','Vtransform','Vstretching','theta_s','theta_b','Tcline','hc','s_rho','s_w','Cs_r','Cs_w','ocean_time'};
  for i=1:numel(var_name)
    ncwrite(nc_ini_new,var_name{i},ncread(nc_ini_old,var_name{i}));
  end
  clear i var_name;
  
  % Transfer grid variables from grid file
  var_name = {{'lon','lat'},{'rho','u','v'}};
  for i=1:numel(var_name{1})
  for j=1:numel(var_name{2})
    var_str = [var_name{1}{i} '_' var_name{2}{j}];  
    ncwrite(nc_ini_new,var_str,ncread(nc_grd_new,var_str));
  end
  end
  clear i j var_name var_str;
  
  % Transfer gridded variables
  var_name = {'zeta','ubar','vbar','u','v','temp','salt','rho'};
  var_dims = [2,2,2,3,3,3,3,3];
  var_grid = {'rho','u','v','u','v','rho','rho','rho'};
  for i=1:numel(var_name)
      
    % Old vars
    x_old = ncread(nc_grd_old,['lon_' var_grid{i}]);
    y_old = ncread(nc_grd_old,['lat_' var_grid{i}]);
    m_old = ncread(nc_grd_old,['mask_', var_grid{i}]);
    z_old = ncread(nc_ini_old,var_name{i});
      
    % Fill value?
    maxval  = max(z_old(abs(z_old)<1e30));
    minval  = min(z_old(abs(z_old)<1e30));
    fillval = unique(z_old(abs(z_old)>1e30));
    if(isempty(fillval)); fillval = 0; end
    
    % New grid
    x_new = ncread(nc_grd_new,['lon_' var_grid{i}]);
    y_new = ncread(nc_grd_new,['lat_' var_grid{i}]);
    m_new = ncread(nc_grd_new,['mask_', var_grid{i}]);
      
    % 2D or 3D?
    if(var_dims(i)==2)
      
      % Do
      z_new   = ones(size(x_new,1),size(x_new,2)).*fillval;
      ntrplnt = scatteredInterpolant( double(x_old(m_old==1)), double(y_old(m_old==1)), z_old(m_old==1), 'linear' );
      z_new(m_new==1) = ntrplnt( x_new(m_new==1), y_new(m_new==1) );
      z_new(z_new < minval) = maxval;
      z_new(z_new > maxval) = minval;
      clear ntrplnt;
        
    elseif(var_dims(i)==3)
        
      % Do
      z_new = ones(size(x_new,1),size(x_new,2),size(z_old,3)).*fillval;
      for j=1:size(z_old,3)
        
        tmp_old = z_old(:,:,j);
        ntrplnt = scatteredInterpolant( double(x_old(m_old==1)), double(y_old(m_old==1)), tmp_old(m_old==1), 'linear' );
        tmp_new = z_new(:,:,j);
        tmp_new(m_new==1) = ntrplnt(x_new(m_new==1), y_new(m_new==1));
        tmp_new(tmp_new < minval) = minval;
        tmp_new(tmp_new > maxval) = maxval;
        z_new(:,:,j) = tmp_new;
        clear tmp_old ntrplnt tmp_new;
          
      end
        
    end
    
    % Save
    ncwrite(nc_ini_new,var_name{i},z_new);
    clear x_old y_old m_old z_old x_new y_new m_new z_new;
      
  end
  clear i;
  
end