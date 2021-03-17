function status = roms_interp_tides(nc_tides_old,nc_tides_new,nc_grd_new)

  % Get info from old file
  info = ncinfo(nc_tides_old);

  % New grid dimensions
  xr_new = ncread(nc_grd_new,'lon_rho');
  xi  = size(xr_new,1);
  eta = size(xr_new,2);
  dim_name  = {info.Dimensions.Name};
  dim_str   = {'xi_rho','eta_rho'};
  dim_lngth = [xi,eta];
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
  if(exist(nc_tides_new,'file')==2); delete(nc_tides_new); end
  ncwriteschema(nc_tides_new,info);
  
  % Write out non-dimensional variables
  for i=[4,11]; ncwrite(nc_tides_new,info.Variables(i).Name,ncread(nc_tides_old,info.Variables(i).Name)); end
  
  % Write out grid variables
  var_name = {'lon_rho','lat_rho','mask_rho'};
  for i=1:numel(var_name)
    ncwrite(nc_tides_new,var_name{i},ncread(nc_grd_new,var_name{i}));
  end

  % Interpolate variables
  xr_old = ncread(nc_tides_old,'lon_rho');
  yr_old = ncread(nc_tides_old,'lat_rho');
  mr_old = ncread(nc_tides_old,'mask_rho');
  xr_new = ncread(nc_tides_new,'lon_rho');
  yr_new = ncread(nc_tides_new,'lat_rho');
  mr_new = ncread(nc_tides_new,'mask_rho');
  var_name = {'tide_Ephase','tide_Eamp','tide_Cphase','tide_Cangle','tide_Cmin','tide_Cmax'};
  for i=1:numel(var_name)
      
    % Old var
    z_old = ncread(nc_tides_old,var_name{i});
    nt    = size(z_old,3);
    
    % Define
    z_new = zeros(size(xr_new,1),size(xr_new,2),nt);
    
    % Loop
    for j=1:nt
      zz = z_old(:,:,j);
      ntrplnt = scatteredInterpolant( double(xr_old(mr_old==1)), double(yr_old(mr_old==1)), double(zz(mr_old==1)), 'linear');
      zn = z_new(:,:,j);
      zn(mr_new==1) = ntrplnt(xr_new(mr_new==1),yr_new(mr_new==1));
    end
      
    % Write
    ncwrite(nc_tides_new,var_name{i},z_new);
    
    % Clean-up
    clear z_old nt z_new zz zn ntrplnt j;
      
  end
  clear i;
  

end