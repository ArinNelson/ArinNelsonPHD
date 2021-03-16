function status = roms_nctruncate(ncname,ix,iy,iz,it)
% status = roms_nctruncate(nc_name,ix,iy)
% truncates gridded variables in nc_name.nc to horizontal indices ix and iy.
% This does take the offests for the u-, v-, and psi-grids into account.
% 
%=========================================================================%
% By Arin Nelson
% On 03/13/2021
%=========================================================================%

  % If fails, status will be 0
  status = 0;

  % Get info
  if(numel(ncname)>3)
  if(strcmp(ncname(end-2:end),'.nc')==1)
      ncname = ncname(1:end-3);
  end
  end
  ncold = [ncname '.nc'];
  info = ncinfo(ncold);
  
  % Parse out some info's
  info_dim_names = {info.Dimensions.Name};  n_dim = numel(info_dim_names);
  info_var_names = {info.Variables.Name};   n_var = numel(info_var_names);
  
  % Search for xi-grid indices
  if(~isempty(ix))
    dim_x_name = {'xi_rho','xi_u','xi_v','xi_psi'};
    dim_x_indx = zeros(4,1);
    for i=1:4
      ii = find( strcmp(dim_x_name{i}, info_dim_names)==1, 1 );
      if(~isempty(ii)); dim_x_indx(i) = ii; end
    end
    do_x = 1;
  else
    do_x = 0;
  end
  
  % Search for eta-grid indices
  if(~isempty(iy))
    dim_y_name = {'eta_rho','eta_u','eta_v','eta_psi'};
    dim_y_indx = zeros(4,1);
    for i=1:4
      ii = find( strcmp(dim_y_name{i}, info_dim_names)==1, 1 );
      if(~isempty(ii)); dim_y_indx(i) = ii; end
    end
    do_y = 1;
  else
    do_y = 0;
  end
  
  % Search for s-grid indices
  if(~isempty(iz))
    dim_z_name = {'s_rho','s_w'};
    dim_z_indx = zeros(2,1);
    for i=1:2
      ii = find( strcmp(dim_z_name{i}, info_dim_names)==1, 1 );
      if(~isempty(ii)); dim_z_indx(i) = ii; end
    end
    do_z = 1;
  else
    do_z = 0;
  end
  
  % Search for a time dimension
  if(~isempty(it))
      
    % Find the (assumedly-single) time dimension  
    dim_t_indx = 0; 
    dim_t_name = '';
    for i=1:n_dim
    if( numel(info.Dimensions(i).Name)>5 )
      if( strcmp(info.Dimensions(i).Name(end-4:end), '_time')==1 )
          dim_t_indx = i; 
          dim_t_name = info.Dimensions(i).Name;
      end
    end
    end
    
    % If no time dimension found, return an error?
    if(dim_t_indx==0 || isempty(dim_t_name))
      error('No time dimension found! Either it doesn''t end in _time or no time variable exists!'); 
    else
      do_t = 1;
    end
    
  else
    do_t = 0;
  end

  % Define new dimension lengths
%   nx = numel(ix);   dim_x_lngth = [nx nx-1 nx nx-1];
%   ny = numel(iy);   dim_y_lngth = [ny ny ny-1 ny-1];
%   nz = numel(iz);   dim_z_lngth = [nz-1 nz];
%   nt = numel(it);   dim_t_lngth = nt;
  n         = do_x*4 + do_y*4 + do_z*2 + do_t;
  dim_name  = cell(n,1);
  dim_lngth = zeros(n,1);
  dim_indx  = zeros(n,1);
  dim_trunc = cell(n,1);
  i         = 1;
  if(do_x) 
    ii = i:(i+3);
    nx = numel(ix);
    dim_name(ii)  = dim_x_name;
    dim_lngth(ii) = [nx nx-1 nx nx-1];
    dim_trunc(ii) = {ix,ix(1:end-1),ix,ix(1:end-1)};
    dim_indx(ii)  = dim_x_indx;
    i  = i + 4;
  end
  if(do_y)
    ii = i:(i+3);
    ny = numel(iy);
    dim_name(ii)  = dim_y_name;
    dim_lngth(ii) = [ny ny ny-1 ny-1];
    dim_trunc(ii) = {iy,iy,iy(1:end-1),iy(1:end-1)};
    dim_indx(ii)  = dim_y_indx;
    i  = i + 4;
  end
  if(do_z)
    ii = i:(i+1);
    nz = numel(iz);
    dim_name(ii)  = dim_z_name;
    dim_lngth(ii) = [nz nz+1];
    dim_trunc(ii) = {iz,[iz(:); iz(end)+1]};
    dim_indx(ii)  = dim_z_indx;
    i  = i + 2;
  end
  if(do_t)
    ii = i;
    nt = numel(it);
    dim_name(ii)  = dim_t_name;
    dim_lngth(ii) = nt;
    dim_trunc{ii} = it;
    dim_indx(ii)  = dim_t_indx;
  end
    
  % Set new dimension lengths in info structure
  for i=1:n
    info.Dimensions(dim_indx(i)).Length = dim_lngth(i);
  end
  
  % It's a little trickier for variable sizes
  for i=1:n_var

    % This variables' dimensions
    var_dims = info.Variables(i).Dimensions;
    
    % Continue if non-monotonic
    if(~isempty(var_dims))
 
      % Find which dimensions are grid variables
      for j=1:numel(var_dims)
        jj = find(strcmp(var_dims(j).Name,{info.Dimensions.Name})==1);
        if(~isempty(jj))
          info.Variables(i).Dimensions(j).Length = info.Dimensions(jj).Length;
        end
        info.Variables(i).Size(j) = info.Variables(i).Dimensions(j).Length;
      end
      clear j jj;

    end 
    clear var_dims;
      
  end
  clear i;
  
  % Create the truncated variable file
  ncnew = [ncname '_trunc.nc'];
  if(exist(ncnew,'file')==2); delete(ncnew); end
  ncwriteschema(ncnew,info);
  
  %-----------------------------------------------------------------------%
  % Transfer variables
  for i=1:numel(info.Variables)
      
    % This variables' dimensions
    var_dims = info.Variables(i).Dimensions;
    
    % If variable is monotonic, transfer it over.
    if(sum(info.Variables(i).Size)==1)
        
      % Transfer the non-gridded variable  
      ncwrite(ncnew,info.Variables(i).Name,ncread(ncold,info.Variables(i).Name)); 
        
    else
        
	  % Get variable size
	  var_size = info.Variables(i).Size;
	  var_indx = squeeze(cell([1,size(var_size)]));
	  var_indx(:) = {':'};
        
      % Check if any dimensions are grid dimensions
      for j=1:numel(var_dims)
        
        % Possible dimensions
        ii = find( strcmp(var_dims(j).Name,dim_name)==1 , 1 );
        if(~isempty(ii))
          var_indx{j} = dim_trunc{ii};
        end
        
      end
      clear j ii;

      % Read in variable
      var_in = ncread(ncold,info.Variables(i).Name);
        
      % Truncate
      var_out = var_in(var_indx{:});
      
      % Write out variable
      ncwrite(ncnew,info.Variables(i).Name,var_out);
        
      % Clean-up
      clear var_in var_out;
      
    end
    clear var_dims;

  end
  clear i;
  
  % If finished, status=1
  status = 1;
  
end