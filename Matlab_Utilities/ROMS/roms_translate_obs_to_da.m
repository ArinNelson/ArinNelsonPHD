function status = roms_translate_obs_to_da(obs_info,mdl_info)
% status = roms_translate_obs_to_da(obs_info,mdl_info)
% Translates the observations with information saved in the structure
% 'obs-info' to the ROMS data assimilation - ready format, which orders the
% observations as a 1D array as follows:
%  1: survey_time
%  2: state_variable
%  3: obs x/y loc
%  4: obs z loc
% 
% obs_info contains the following info:
% obsFile  : cell array of NetCDF files containing obs from each inst.
% obsLabel : cell array of names/labels relating to each instrument
% obsXI    : vector of fractional x-positions on model grid
%          : Note: grid cell centers are xi=0.5:1:L-0.5
% obsYI    : vector of fractional y-positinos on model grid
%            Note: grid cell centers are yi=0.5:1:M-0.5
% obsDepth : cell array of depths measured by each instrument
%            Note: Negative values = actual depth
%                  Positive values = ROMS depth level
% obsTime  : cell array of times sampled by each instrument
%            Note: must be in units of days since reference time in roms.in
% obsVars  : cell array of strings listing variables present in each file
% 
% Optional contents
% obsLon   : longitudinal positions of each instrument
% obsLat   : latitudinal positions of each instrument
% obsVarTranslate : (optinal) relates variable names from obs file(s) to
%                   ROMS variable names

%==========================================================================
% by Arin Nelson
% on 02/10/2021
%==========================================================================

  % If it makes it to the end, status=1
  status = 1;

  % Validate inputs
  if(~iscell(obs_info.obsTime));    obs_info.obsTime = {obs_info.obsTime};  end
  
  % Initialize the output file's schema
  schema = roms_netcdf_template_obs;

  % Set attributes as defined by the obs_info
  schema.atts{2,2} = mdl_info.Title;
  schema.atts{4,2} = mdl_info.Label;        %n_obs = numel(mdl_info.Label);
  schema.atts{5,2} = 'zeta,ubar,vbar,u,v,temp,salt';
  schema.atts{6,2} = 'riroms_grd_thirded.nc';
  schema.atts{10,2} = mdl_info.LM;

  % Obs names into a single string
  schema.atts{8,2} = ''; 
  for i=1:(numel(obs_info.obsLabel)-1)
    schema.atts{8,2} = [schema.atts{8,2} obs_info.obsLabel{i} ', '];
  end
  schema.atts{8,2} = [schema.atts{8,2} obs_info.obsLabel{end}];
   
  % Determine # surveys
  t = [];
  for i=1:numel(obs_info.obsTime)
    t = [t(:); obs_info.obsTime{i}(:)];
  end
  survey_time = unique(t);
  n_survey    = numel(survey_time);
  schema.dims{1,2} = n_survey;
   
%   % Determine # state vars
%   if(~iscell(obs_info.obsVars))
%     var_names = strsplit(obs_info.obsVars,{' ',','});
%   else
%     var_names = obs_info.obsVars;    
%   end
%   n_var    = numel(var_names);
  schema.dims{2,2} = 7;
  
  % Determine # data points
  n_datum = 5*7*n_survey;   % SIMPLE FOR NOW
  schema.dims{3,2} = n_datum;
  
  % With all dimensions known, the NetCDF file can be initialized
  fOut = 'obs_combo.nc';        % SIMPLE FOR NOW
  roms_ncgen(fOut,schema);
  
  % Determine indices for data to gather
  % ROMS documentation recommends this order:
  % 1: survey_time
  % 2: state_variable
  % 3: obs x/y loc
  % 4: obs z loc
  % (FOR NOW, MAKING THIS SIMPLE SO I CAN GET WORKING WITH IT!)
  n = 1;
  for is=1:n_survey
    for j=1:5  
    for io=1:7
          
        % Get data (SIMPLE FOR NOW)  
            if(j==1)   
                dat = ncread(obs_info.obsFile{io},'zeta',is,1); 
                ot  = 1;
                vz  = 1;
            end
            if(j==2) 
                dat = ncread(obs_info.obsFile{io},'temp',[1 is],[1 1]); 
                ot  = 6;
                vz  = 1;
            end
            if(j==3)   
                dat = ncread(obs_info.obsFile{io},'temp',[2 is],[1 1]); 
                ot  = 6;
                vz  = 2;
            end
            if(j==4)   
                dat = ncread(obs_info.obsFile{io},'salt',[1 is],[1 1]); 
                ot  = 7;
                vz  = 1;
            end
            if(j==5)   
                dat = ncread(obs_info.obsFile{io},'salt',[2 is],[1 1]); 
                ot  = 7;
                vz  = 2;
            end
        
          % Write to file (SIMPLE FOR NOW
          ncwrite(fOut,'obs_type',ot,n);
          ncwrite(fOut,'obs_provenance',io,n);
          ncwrite(fOut,'obs_time',survey_time(is),n);
          ncwrite(fOut,'obs_depth',obs_info.obsZI(vz),n);
          ncwrite(fOut,'obs_Xgrid',obs_info.obsXI(io)*(330/1000),n);
          ncwrite(fOut,'obs_Ygrid',obs_info.obsYI(io)*(360/1100),n);
          ncwrite(fOut,'obs_Zgrid',obs_info.obsZI(vz),n);
          ncwrite(fOut,'obs_error',0.001,n);
          ncwrite(fOut,'obs_value',dat,n);
          n = n + 1;
       
      end
      end
  end
  clear is io;
 
  % Done!
  
end