function status = roms_sample_output(mdl_info, obs_info)
% status = roms_sample_output(mdl_info, obs_info, output_file)
% Take synthetic observations from the simulation with info in structure
% 'mdl_info' based on information from the observation platforms in
% structure 'obs_info'.  The synthetic observations are saved into the
% NetCDF file 'output_file'.
% 
%==========================================================================
% by Arin Nelson
% on 02/08/2021
%==========================================================================

  % If it makes it to the end, status = 1
  status = 1;

  % Check & Validate Inputs
  if(~iscell(obs_info.Label));  obs_info.Label = {obs_info.Label};  end
  if(~iscell(obs_info.Time));   obs_info.Time  = {obs_info.Time};   end
  if(~iscell(obs_info.Depth));  obs_info.Depth = {obs_info.Depth};  end
  if(~iscell(obs_info.ZI));     obs_info.ZI    = {obs_info.ZI};     end
  if(~iscell(obs_info.TI));     obs_info.TI    = {obs_info.TI};     end
  n_obs_locs = numel(obs_info.Label);

  % Options that may be changed later
  tol = 1e5;        % Tolerance of defining an index as an integer value
  
  % List model files
  mdl_filelist = ls([mdl_info.Dir '/' mdl_info.Prfx '*.nc']);

  % Get number of time steps per file
  nt_per_file = numel(ncread([mdl_info.Dir '/' mdl_filelist(1,:)],'ocean_time'));

%   % Get time information
%   n_files     = size(mdl_filelist,1);
%   nt          = n_files*nt_per_file;
%   mdl_time_value     = NaN(nt,1);
%   mdl_time_fileindex = NaN(nt,1);

%   % Read in ocean_time from each file
%   n = 1;
%   for i=1:n_files
%     ii = n:(n+nt_per_file-1);  
%     mdl_time_value(ii)     = ncread([mdl_info.Dir '/' mdl_filelist(i,:)],'ocean_time');
%     mdl_time_fileindex(ii) = i;
%     n = n + nt_per_file;
%   end
%   clear nt n i ii;
  
%   % Transform depends on file... switch to using indices ifile and iinfile
%   mdl_time_value = mdl_time_value./(5760*15) + datenum('2016-01-01');
  
%   % Determine which files will be needed for the synthetic observations
%   time_range = [];
%   if(~iscell(obs_info.Time))
%     time_range = obs_info.Time([1 end]);
%   else
%     time_range = obs_info.Time{1}([1 end]);
%     for i=2:numel(obs_info.Time)
%       time_range(1) = min([time_range(1), obs_info.Time{i}(1)]);
%       time_range(2) = max([time_range(2), obs_info.Time{i}(2)]);
%     end
%   end
%   
%   % Only keep list of files available within this timeframe
%   i1 = find(mdl_time_value<=time_range(1),1,'last');    if(isempty(i1));    i1=1;           end
%   i2 = find(mdl_time_value>=time_range(2),1,'first');   if(isempty(i2));	i2=n_files;     end
%   ii = i1:i2;
%   mdl_time_value = mdl_time_value(ii);
%   mdl_time_fileindex = mdl_time_fileindex(ii);
%   mdl_filelist(~ismember(1:n_files,unique(mdl_time_fileindex)),:) = [];
%   mdl_time_fileindex = mdl_time_fileindex - min(mdl_time_fileindex)+1;
  
%   % Read in model grid information
%   [xr,yr,mr] = roms_read_grid(mdl_info.Grid,'rho','lon');
%   [xu,yu,mu] = roms_read_grid(mdl_info.Grid,'u',  'lon');
%   [xv,yv,mv] = roms_read_grid(mdl_info.Grid,'v',  'lon');
%   [xp,yp,mp] = roms_read_grid(mdl_info.Grid,'psi','lon');
%   mdl_grid_types = {'rho','u','v','psi'};
%   mdl_lon  = {xr,xu,xv,xp};     clear xr xu xv xp;
%   mdl_lat  = {yr,yu,yv,yp};     clear yr yu yv yp;
%   mdl_mask = {mr,mu,mv,mp};     clear mr mu mv mp;
  
  % Since we already computed the x,y,z,t indices, we can perform the
  % sampling loop in almost any order with equal efficiency.  I will go
  % through in order as follows:
  % Observation Station - Observed Variable - within XYZT range
  for i=1:n_obs_locs
      
    % Get location and time info from this observing station
    xi = obs_info.XI(i);
    yi = obs_info.YI(i);
    %t  = obs_info.Time{i};  nt = numel(t);
    if(numel(obs_info.ZI)>1)
      zi = obs_info.ZI{i};  
    else
      zi = obs_info.ZI{1};     
    end
    nz = numel(zi);
    if(numel(obs_info.TI)>1)
      ti = obs_info.TI{i};    
      tt = obs_info.Time{i};
    else
      ti = obs_info.TI{1};
      tt = obs_info.Time{1};
    end
    nt = numel(ti);
    
    % Get list of variables that will be read in  
    if(n_obs_locs == 1)
      if(~iscell(obs_info.Vars));   v = strsplit(obs_info.Vars,{' ',','});
      else;                         v = obs_info.Vars;
      end
    else
      if(~iscell(obs_info.Vars));   v = strsplit(obs_info.Vars,{' ',','});
      else;                         v = strsplit(obs_info.Vars{i},{' ',','});
      end
    end
    
    % Determine if each variable is 2D or 3D in the first output file
    info = ncinfo([mdl_info.Dir '/' mdl_filelist(1,:)]);
    d    = zeros(numel(v),1);  % # dimensions of variable (2D or 3D)
    for j=1:numel(v)
      d(j) = numel(info.Variables( strcmp({info.Variables.Name},v{j})==1 ).Dimensions)-1;  % Minus the time dimension
    end
    clear info j;
    
    % Construct the observation file name
    fOut = [mdl_info.Label '_obs_' obs_info.Label{i} '.nc'];
    
    % Generate the observation file
    ncid = netcdf.create(fOut,'CLOBBER');
    dimz = netcdf.defDim(ncid,'nz',nz);
    dimt = netcdf.defDim(ncid,'nt',nt);
    netcdf.defVar(ncid,'lon','NC_DOUBLE',[]);
    netcdf.defVar(ncid,'lat','NC_DOUBLE',[]);
    netcdf.defVar(ncid,'time','NC_DOUBLE',dimt);
    netcdf.defVar(ncid,'x_index','NC_DOUBLE',[]);
    netcdf.defVar(ncid,'y_index','NC_DOUBLE',[]);
    netcdf.defVar(ncid,'z_index','NC_DOUBLE',dimz);
    netcdf.defVar(ncid,'t_index','NC_DOUBLE',dimt);
    for j=1:numel(v)
        
      % Determine if given variable is 2D or 3D in output
      % (unless there's only 1 z level sampled, then it doesn't matter...)  
      if(nz==1)
        netcdf.defVar(ncid,v{j},'NC_DOUBLE',[dimt]);
      else
        if(d(j)==2)
          netcdf.defVar(ncid,v{j},'NC_DOUBLE',[dimt]);
        elseif(d(j)==3)
          netcdf.defVar(ncid,v{j},'NC_DOUBLE',[dimz dimt]);  
        else
          error(['Variable ' v{j} ' has an unsupported number of spatial dimensions: ' num2str(d(j))]);
        end
      end
      
    end
    netcdf.endDef(ncid);
    netcdf.close(ncid);
    clear j dimz dimt ncid;
    
    % Write some initially-known things
    ncwrite(fOut,'lon',obs_info.Lon(i));
    ncwrite(fOut,'lat',obs_info.Lat(i));
    ncwrite(fOut,'time',tt);
    ncwrite(fOut,'x_index',xi);
    ncwrite(fOut,'y_index',yi);
    ncwrite(fOut,'z_index',zi);
    ncwrite(fOut,'t_index',ti);
    
    % Loop through variables
    for j=1:numel(v)
        
      % Gather indices around points of interest
      if( abs(xi-round(xi)) < (1/tol) );	ix = round(xi);
      else;                                 ix = [floor(xi) ceil(xi)];
      end
      
      if( abs(yi-round(yi)) < (1/tol) );    iy = round(yi);
      else;                                 iy = [floor(yi) ceil(yi)];
      end
      
      if(d(j)==2)
        iz = 1;
      elseif(nz==1)
        if( abs(zi-round(zi)) < (1/tol) );  iz = round(zi);
        else;                               iz = [floor(zi) ceil(zi)];
        end
      else
        iz = [];
        for k=1:nz
          if( abs(zi(k)-round(zi(k))) < (1/tol) )
            iz = [iz(:); round(zi(k))];
          else
            iz = [iz(:); floor(zi(k)); ceil(zi(k))];
          end
        end
      end
      iz = unique(iz(~isnan(iz)));
      
      it = [];
      for k=1:nt
        if( abs(ti(k)-round(ti(k))) < (1/tol) )
          it = [it(:); round(ti(k))];
        else
          it = [it(:); floor(ti(k)); ceil(ti(k))];
        end
      end
      it = unique(it(~isnan(it)));
      
      % Load variable at the given indices for all time steps
      ifile   = ceil(it./nt_per_file);
      iinfile = mod(it,nt_per_file);    iinfile(iinfile==0) = nt_per_file;
      tmp_var = zeros(numel(ix),numel(iy),numel(iz),numel(it));
      for k=1:numel(ifile)
      if(~isnan(ifile(k)))  
      for m=1:numel(iz)
      if(d(j)==2)    
        tmp_var(:,:,m,k) = ncread([mdl_info.Dir '/' mdl_filelist(ifile(k),:)],v{j},[ix(1),iy(1),iinfile(k)],[numel(ix),numel(iy),1]);
      else
        tmp_var(:,:,m,k) = ncread([mdl_info.Dir '/' mdl_filelist(ifile(k),:)],v{j},[ix(1),iy(1),iz(m),iinfile(k)],[numel(ix),numel(iy),1,1]);
      end 
      end
      end
      end
      
      % Perform interpolations in space and time
      ii = {ix,iy,iz,it};
      xx = {xi,yi,zi,ti};
      ni = zeros(4,1);  for k=1:4;  ni(k)=numel(ii{k}); end
      ii(ni==1) = [];
      xx(ni==1) = [];
      zz = squeeze(tmp_var);
      
      % Depending on the number of non-1-length dimensions determines which
      % interpolation method to use
      if(numel(ii)==4)
        ntrpd_var = interpn(ix,iy,iz,it,zz,xi,yi,zi,ti,'linear');
      elseif(numel(ii)==3)
        ntrpd_var = interpn(ii{1},ii{2},ii{3},zz,xx{1},xx{2},xx{3},'linear');
      elseif(numel(ii)==2)
        ntrpd_var = interp2(ii{1},ii{2},zz,xx{1},xx{2},'linear');
      elseif(numel(ii)==1)
        ntrpd_var = interp1(ii{1},zz,xx{1},'linear');
      end
      ntrpd_var = squeeze(ntrpd_var);
        
      % Write to file
      ncwrite(fOut,v{j},ntrpd_var);
      
      % Clean-up
      clear ix iy iz it ifile iinfile tmp_var ntrpd_var;
      
    end
    clear j;
    
    % Clean-up
    clear xi yi zi ti d;
      
  end
  clear i;
  
  
end