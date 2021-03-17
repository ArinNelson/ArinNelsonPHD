function status = roms_interp_bry(nc_bry_old,nc_grd_old,nc_bry_new,nc_grd_new)

  % Get info from old file
  info = ncinfo(nc_bry_old);

  % New grid dimensions
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
  if(exist(nc_bry_new,'file')==2); delete(nc_bry_new); end
  ncwriteschema(nc_bry_new,info);
  
  % Write out non-dimensional variables
  for i=[1:11,30]; ncwrite(nc_bry_new,info.Variables(i).Name,ncread(nc_bry_old,info.Variables(i).Name)); end
  
  % Write out new grid variables
  dir_str = {'west','south','east','north'};
  var_str = {'lon','lat'};
  var_grd = {'rho','u','v','psi'};
  for i=1:numel(var_str)
  for j=1:numel(var_grd)
  for k=1:numel(dir_str)
      
    % Var
    var_name = [var_str{i} '_' var_grd{j} '_' dir_str{k}];
    ijk = find( strcmp({info.Variables.Name},var_name)==1 );
    if(~isempty(ijk))
        
      % Load 'n' save
      var_tmp = ncread(nc_grd_new,[var_str{i} '_' var_grd{j}]);
      switch k
          case 1;   var_tmp = var_tmp(1,:);
          case 2;   var_tmp = var_tmp(:,1);
          case 3;   var_tmp = var_tmp(end,:);
          case 4;   var_tmp = var_tmp(:,end);
      end
      ncwrite(nc_bry_new,var_name,var_tmp);  
      
    end
    
  end
  end
  end
  clear i j k;
  
  % Interpolate variables
  var_name = {'zeta','ubar','vbar','u','v','temp','salt'};
  var_dims = [2 2 2 3 3 3 3];
  var_dir  = {'west','south','east','north'};
  var_dim  = {'lat','lon','lat','lon'};
  var_grd = {'rho','u','v','u','v','rho','rho'};
  for i=1:numel(var_name)
  for j=1:numel(var_dir)
      
    % Continue if variable exists
    ij = find( strcmp([var_name{i} '_' var_dir{j}],{info.Variables.Name})==1 );
    if(~isempty(ij))
      
      % Grid variables
      x_old = ncread(nc_bry_old,[var_dim{j} '_' var_grd{i} '_' var_dir{j}]);
      z_old = ncread(nc_bry_old,[var_name{i} '_' var_dir{j}]);
      x_new = ncread(nc_bry_new,[var_dim{j} '_' var_grd{i} '_' var_dir{j}]); 
      
      % And masks
      m_old = ncread(nc_grd_old,['mask_' var_grd{i}]);
      m_new = ncread(nc_grd_new,['mask_' var_grd{i}]);
      switch j
          case 1;   m_old = m_old(1,:);     m_new = m_new(1,:);
          case 2;   m_old = m_old(:,1);     m_new = m_new(:,1);
          case 3;   m_old = m_old(end,:);   m_new = m_new(end,:);
          case 4;   m_old = m_old(:,end);   m_new = m_new(:,end);
      end
      
      % Interp
      if(var_dims(i)==2)
        z_new = zeros(numel(x_new),size(z_old,2));
        z_new(m_new==1,:) = interp1(x_old(m_old==1),z_old(m_old==1,:),x_new(m_new==1),'linear');
      elseif(var_dims(i)==3)
        z_new = zeros(numel(x_new),size(z_old,2),size(z_old,3));
        z_new(m_new==1,:,:) = interp1(x_old(m_old==1),z_old(m_old==1,:,:),x_new(m_new==1),'linear');
      end
      
      % Write
      ncwrite(nc_bry_new,[var_name{i} '_' var_dir{j}],z_new);
      
      % Clean-up
      clear x_old z_old x_new m_old m_new z_new;
        
    end
      
  end
  end
  clear i j;

end