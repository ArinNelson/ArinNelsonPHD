function status = roms_sample_output_roughversion(output_dir, output_prefix, obs_info)
%=========================================================================%
% status = roms_sample_output(output_dir, output_prefix, var_list, obs_info)
% Sample the ROMS output in directory 'output_dir' using the observation
% information in structure 'obs_info'.  The output files are expected to
% have the format 'output_prefix'_####.nc, where #### is the output
% timestamp.
% 
% The structure 'obs_info' contains the following information:
%  Title         : string to put as 'title' attribute in the output file
%  FileName      : save observations into 'FileName'.nc            (string)
%  GridDims      : [x,y,z] dimensions of grid                  (1x3 vector)
%  GridFile      : grid file associated with this ROMS output      (string)
%  GridSpherical : spherical switch w/in grid file                 (string)
%  ObsName       : list of instruments/labels/etc. of observation locations
%  ObsVars       : Names of variables sampled by each instrument (cell arra
%  ObsT          : Time(s) of instrument samples    (cell array of vectors)
%  ObsX          : x-position of instrument(s) (either m or deg.east)      
%  ObsXI         : fractional x-index on grid
%  ObsY          : y-position of instrument(s) (either m or deg.north)
%  ObsYI         : fractional y-index on grid
%  ObsZ          : Depth(s) of instrument samples   (cell array of vectors)
%  StateVars     : list of unique state var.s that samples will be taken of
%  XType         : horizontal grid type (either 'lon' or 'x')      (string)
% 
%=========================================================================%
% By Arin Nelson on 02/08/2021
%=========================================================================%

  % Validate inputs

  
  % If it makes it to the end, status=1
  status = 1;
  
%   % Read in info from grid file
%   [x_rho,y_rho] = roms_read_grid(obs_info.GridFile,'rho',obs_info.XType);
%   [x_u,y_u]     = roms_read_grid(obs_info.GridFile,'u',  obs_info.XType);
%   [x_v,y_v]     = roms_read_grid(obs_info.GridFile,'v',  obs_info.XType);
%   [x_psi,y_psi] = roms_read_grid(obs_info.GridFile,'psi',obs_info.XType);
%   x = {x_rho,x_u,x_v,x_psi};
%   y = {y_rho,y_u,y_v,y_psi};
  
%   % Constant: variables' grid type
%   gridtype_str   = {'zeta','ubar','vbar','u','v','temp','salt'};
%   gridtype_index = [1 2 3 2 3 1 1]; %rho, u, v
  
  % Initialize the output file's schema
  schema = roms_netcdf_template_obs;
  
  % Set attributes as defined by the obs_info
  schema.atts{2,2} = obs_info.Title;
  schema.atts{4,2} = obs_info.GridFile;
  schema.atts{5,2} = obs_info.StateVars;
  schema.atts{6,2} = output_dir;
  schema.atts{10,2} = obs_info.GridDims;
  
  % Obs names into a single string
   schema.atts{8,2} = ''; 
   for i=1:(numel(obs_info.ObsName)-1)
     schema.atts{8,2} = [schema.atts{8,2} obs_info.ObsName{i} ', '];
   end
   schema.atts{8,2} = [schema.atts{8,2} obs_info.ObsName{end}];
  
  % Determine # surveys
  n_obs = numel(obs_info.ObsName);
  t = [];
  for i=1:n_obs
    t = [t(:); obs_info.ObsT{i}(:)];
  end
  survey_time = unique(t);
  n_survey    = numel(survey_time);
  schema.dims{1,2} = n_survey;
  
  % Determine # state vars
  var_names = strsplit(obs_info.StateVars,{' ',','});
  n_var    = numel(var_names);
  schema.dims{2,2} = n_var;
  
  % Determine # datums per instrument
  n_inst = numel(obs_info.ObsName);
  inst_x = zeros(n_inst,1);
  inst_y = zeros(n_inst,1);
  inst_z = cell(n_inst,1);
  inst_t = cell(n_inst,1);
  inst_n = zeros(n_inst,1);
  inst_v = cell(n_inst,1);
  for i=1:n_inst
      
    % Get instruments' sampling info 
    inst_x(i) = obs_info.ObsX(i);
    inst_y(i) = obs_info.ObsY(i);
    inst_z{i} = obs_info.ObsZI{i};
    inst_t{i} = obs_info.ObsT{i};
    inst_v{i} = obs_info.ObsVars{i};
    
    % Total number of data points per instrument is 
    %       # depths sampled 
    % times # times sampled 
    % times # variables sampled
    inst_n(i) = numel(inst_z{i}) * numel(inst_t{i}) * numel(inst_v{i});    
    
  end
  clear i;
  
  % Number of observations per survey time
  n_obs = zeros(n_survey,1);
  for i=1:n_inst
    tmp = inst_t{i};
    nz  = numel(inst_z{i});
    nv  = numel(inst_v{i});
    for j=1:numel(tmp)
      jj = find(survey_time == tmp(j));  
      n_obs(jj) = n_obs(jj) + nz*nv;
    end
  end
  clear i j tmp nz jj;

  % The total number of datums
  n_datum = sum(n_obs);
  schema.dims{3,2} = n_datum;
  
  % With all dimensions known, the NetCDF file can be initialized
  fOut = [obs_info.FileName '.nc'];
  roms_ncgen(fOut,schema);
  
  % Write some variables
  ncwrite(fOut,'spherical',obs_info.GridSpherical);
  ncwrite(fOut,'Nobs',n_datum);
  ncwrite(fOut,'survey_time',survey_time);
  
  % Sim times
  sim_time = roms_read_output_times(output_dir,output_prefix);
  
  % Generate the index array corresponding to the observations
  % 1: survey_time
  % 2: state_variable
  % 3: obs x/y loc
  % 4: obs z loc
  n         = 1;
  %obs_index = zeros(n_datum,4);
  for i=1:n_survey
  isamptime = find(sim_time==survey_time(i)); 
      
    % Loop through obs locs
    jj = [];
    for j=1:numel(obs_info.ObsName)
      if(any( obs_info.ObsT{j} == survey_time(i) )) 
        jj = [jj(:), j];
      end
    end
    
    % Get variables and depth levels read in by each instrument at this time
    if(~isempty(jj))
      obs_var_info = cell(numel(jj),1);
      for j=1:numel(jj)
        tmp = ''; for k=1:numel(obs_info.ObsVars{jj(j)}); tmp = [tmp ',' obs_info.ObsVars{jj(j)}{k}]; end
        obs_var_info{j} = [obs_info.ObsName{jj(j)} tmp ',' num2str([obs_info.ObsZI{jj(j)}])];
      end
    end
    
    % Get unique variable and depth combo's
    if(numel(obs_var_info)>1)
      obs_var_unique = unique(obs_var_info);
    else
      obs_var_unique = obs_var_info;
    end
    
    % Gather unique variable and depth level combo's
    unique_vars = cell(0,1);
    for j=1:numel(obs_var_unique)
      tmp = strsplit(obs_var_unique{j},',');
      tmp1 = tmp(2:end-1);
      tmp2 = strsplit(tmp{end},' ');
      for k=1:numel(tmp1)
      for m=1:numel(tmp2)
        unique_vars{end+1} = [tmp1{k} ',' tmp2{m}];
      end
      end
    end
    unique_all = unique(unique_vars);
    clear tmp tmp1 tmp2 j k m;
    
    % Loop through unique variable & depth combo's
    for j=1:numel(unique_all)
        
      % Split into variable name and depth
      tmp = strsplit(unique_all{j},',');
      var_name = tmp{1};
      var_dpth = str2num(tmp{2});
      clear tmp;
        
%       % Get grid type of observation
%       xx = x{gridtype_index( strcmp(gridtype_str,var_name)==1 )};
%       yy = y{gridtype_index( strcmp(gridtype_str,var_name)==1 )};
 
      % Read in this variable at this depth
      zz = roms_read_outputs_at_time(output_dir,output_prefix,{var_name},var_dpth,isamptime);
      %zz = zz{1};
      
      % Loop through instruments that observe this variable at this depth and time
      kk = [];
      if(numel(jj)==1)
        kk = jj;
      else
        for k=1:numel(jj)
          test_depths = obs_info.ObsZ{jj(k)};
          test_vars   = obs_info.ObsVars{jj(k)};
          if(any(test_depths == var_dpth) & any(strcmp(test_vars,var_name)))
            kk = [kk(:), jj(k)];
          end
        end
      end
      clear k test_*
      
      % Output loop
      for k=1:numel(kk)

        % Write to output file
        ncwrite(fOut,'obs_type',find( strcmp(var_names,var_name)==1 ),n);
        ncwrite(fOut,'obs_provenance',kk(k),n);
        ncwrite(fOut,'obs_time',survey_time(i),n);
        ncwrite(fOut,'obs_lon',obs_info.ObsX(kk(k)),n);
        ncwrite(fOut,'obs_lat',obs_info.ObsY(kk(k)),n);
        ncwrite(fOut,'obs_depth',obs_info.ObsZ{kk(k)}(var_dpth),n);
        ncwrite(fOut,'obs_Xgrid',obs_info.ObsXI(kk(k)),n);
        ncwrite(fOut,'obs_Ygrid',obs_info.ObsYI(kk(k)),n);
        ncwrite(fOut,'obs_Zgrid',obs_info.ObsZI{kk(k)}(var_dpth),n);
        ncwrite(fOut,'obs_error',0.001,n);
        ncwrite(fOut,'obs_value',zz(obs_info.ObsXI(kk(k)),obs_info.ObsYI(kk(k))),n);
        
        % Next time step
        n = n + 1;
        
      end
      clear k;
      
      % Clean-up
      clear zz var_name var_dpth 
      
    end
    clear j jj;

  end
  clear i;

end