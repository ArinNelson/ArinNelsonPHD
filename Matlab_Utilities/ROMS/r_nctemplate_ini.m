function info = r_nctemplate_ini(spherical)
%=========================================================================%
% info = r_nctemplate_ini(spherical)
% 
% Generate a structure containing dimension and variable information of a
% typical ROMS initial conditions NetCDF file.  Dimension lengths and
% variable values are initially blank, but can be set manually.
% 
% Alternately, calling r_ncgen_ini will generate a template initial 
% conditions NetCDF file with dimension lengths and variable values present 
% in the roms.in file specified by this function.
%=========================================================================%
% by Arin Nelson on 09/23/2020
%=========================================================================%

% Define the dimension names and indices
info.dims = struct;
info.dims(01).name = 'xi_rho';     
info.dims(02).name = 'xi_u';       
info.dims(03).name = 'xi_v';       
info.dims(04).name = 'xi_psi';     
info.dims(05).name = 'eta_rho';    
info.dims(06).name = 'eta_u';     
info.dims(07).name = 'eta_v';      
info.dims(08).name = 'eta_psi';    
info.dims(09).name = 's_rho';      
info.dims(10).name = 's_w';        
info.dims(11).name = 'tracer';     
info.dims(12).name = 'ocean_time';

% Define the global attributes
info.atts = struct;
info.atts(1).name = 'type';
info.atts(2).name = 'history';

% Set the global attributes
info.atts(1).value = 'INITIALIZATION file';
info.atts(2).value = ['Initialized by r_nctemplate_ini at ' datestr(now) '.'];

% Define the vertical grid configuration variables
info.vars = [];
info.vars(end+1).name = 'spherical';
info.vars(end  ).info = {'type','NC_INT'; ...
                           'dimid',[]; ...
                           'long_name','grid type logical switch'; ...
                           'flag_values',[0 1]; ...
                           'flag_meanings','Cartesian, Spherical'};
info.vars(end+1).name = 'Vtransform';
info.vars(end  ).info = {'type','NC_INT'; ...
                           'dimid',[]; ...
                           'long_name','vertical terrain-following transformation equation'};
info.vars(end+1).name = 'Vstretching';
info.vars(end  ).info = {'type','NC_INT'; ...
                           'dimid',[]; ...
                           'long_name','vertical terrain-following stretching function'};
info.vars(end+1).name = 'theta_s';
info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                           'dimid',[]; ...
                           'long_name','S-coordinate surface control parameter'};                
info.vars(end+1).name = 'theta_b';
info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                           'dimid',[]; ...
                           'long_name','S-coordinate bottom control parameter'};                       
info.vars(end+1).name = 'Tcline';
info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                           'dimid',[]; ...
                           'long_name','S-coordinate surface/bottom layer width'; ...
                           'units','meters'};     
info.vars(end+1).name = 'hc';
info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                           'dimid',[]; ...
                           'long_name','S-coordinate parameter, critical depth'; ...
                           'units','meters'};   
info.vars(end+1).name = 's_rho';
info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                           'dimid',[8]; ...
                           'long_name','S-coordinate at RHO-points'; ...
                           'valid_min',-1; ...
                           'valid_max',0; ...
                           'positive','up'; ...
                           'standard_name','ocean_s_coordinate_g(Vtransform)'; ...
                           'formula_terms','s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc' };
info.vars(end+1).name = 's_w';
info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                           'dimid',[9]; ...
                           'long_name','S-coordinate at W-points'; ...
                           'valid_min',-1; ...
                           'valid_max',0; ...
                           'positive','up'; ...
                           'standard_name','ocean_s_coordinate_g(Vtransform)'; ...
                           'formula_terms','s: s_w C: Cs_w eta: zeta depth: h depth_c: hc' };
info.vars(end+1).name = 'Cs_r';
info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                           'dimid',[8]; ...
                           'long_name','S-coordinate stretching function at RHO-points'; ...
                           'valid_min',-1; ...
                           'valid_max',0 };
info.vars(end+1).name = 'Cs_w';
info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                           'dimid',[9]; ...
                           'long_name','S-coordinate stretching function at W-points'; ...
                           'valid_min',-1; ...
                           'valid_max',0 };
                   
% Define horizontal grid variables independent of spherical switch
info.vars(end+1).name = 'h';
info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                           'dimid',[0 4]; ...
                           'long_name','bathymetry at RHO-points'; ...
                           'units','meters' };
                                   
% Define horizontal grid variables dependent of spherical switch                   
if(spherical==1)   
    
  % Last attribute of h  
  info.vars(end).info(end+1,:) = {'coordinates','lon_rho,lat_rho'};
  
  % Grid variables
  info.vars(end+1).name = 'lon_rho';
  info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                              'dimid',[0 4]; ...
                              'long_name','longitude at RHO-points'; ...
                              'units','degrees_east'; ...
                              'standard_name','longitude' };
  info.vars(end+1).name = 'lat_rho';
  info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                             'dimid',[0 4]; ...
                             'long_name','latitude at RHO-points'; ...
                             'units','degrees_north'; ...
                             'standard_name','latitude' };          
  info.vars(end+1).name = 'lon_u';
  info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                             'dimid',[1 5]; ...
                             'long_name','longitude at U-points'; ...
                             'units','degrees_east'; ...
                             'standard_name','longitude' };
  info.vars(end+1).name = 'lat_u';
  info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                             'dimid',[1 5]; ...
                             'long_name','latitude at U-points'; ...
                             'units','degrees_north'; ...
                             'standard_name','latitude' };                        
  info.vars(end+1).name = 'lon_v';
  info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                             'dimid',[2 6]; ...
                             'long_name','longitude at V-points'; ...
                             'units','degrees_east'; ...
                             'standard_name','longitude' };
  info.vars(end+1).name = 'lat_v';
  info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                             'dimid',[2 6]; ...
                             'long_name','latitude at V-points'; ...
                             'units','degrees_north'; ...
                             'standard_name','latitude' };     
elseif(spherical==0)
    
  % Last info of h  
  info.vars(end).info(end+1,:) = {'coordinates','x_rho,y_rho'};
  
  % Grid variables
  info.vars(end+1).name = 'x_rho';
  info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                             'dimid',[0 4]; ...
                             'long_name','X-location at RHO-points'; ...
                             'units','meters' };
  info.vars(end+1).name = 'y_rho';
  info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                             'dimid',[0 4]; ...
                             'long_name','Y-location at RHO-points'; ...
                             'units','meters' };          
  info.vars(end+1).name = 'x_u';
  info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                             'dimid',[1 5]; ...
                             'long_name','X-location at U-points'; ...
                             'units','meters' };
  info.vars(end+1).name = 'y_u';
  info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                             'dimid',[1 5]; ...
                             'long_name','Y-location at U-points'; ...
                             'units','meters' };                        
  info.vars(end+1).name = 'x_v';
  info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                             'dimid',[2 6]; ...
                             'long_name','X-location at V-points'; ...
                             'units','meters' };
  info.vars(end+1).name = 'y_v';
  info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                             'dimid',[2 6]; ...
                             'long_name','Y-location at V-points'; ...
                             'units','meters' };     
else
  error(['Unknown value for spherical flag: ' num2str(spherical) '. Must be 0 or 1.']);
end

% Define initial condition variables
info.vars(end+1).name = 'ocean_time';
info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                           'dimid',[11]; ...
                           'long_name','time since initialization'; ...
                           'units','seconds' };
info.vars(end+1).name = 'zeta';
info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                           'dimid',[0 4 11]; ...
                           'long_name','free surface'; ...
                           'units','meters'; ...
                           'time','ocean_time' };
info.vars(end+1).name = 'ubar';
info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                           'dimid',[1 5 11]; ...
                           'long_name','vertically integrated u-momentum component'; ...
                           'units','meters seconds-1'; ...
                           'time','ocean_time' };
info.vars(end+1).name = 'vbar';                      
info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                           'dimid',[2 6 11]; ...
                           'long_name','vertically integrated v-momentum component'; ...
                           'units','meters seconds-1'; ...
                           'time','ocean_time' };            
info.vars(end+1).name = 'u';
info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                           'dimid',[1 5 8 11]; ...
                           'long_name','u-momentum component'; ...
                           'units','meters seconds-1'; ...
                           'time','ocean_time' };
info.vars(end+1).name = 'v';
info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                           'dimid',[2 6 8 11]; ...
                           'long_name','v-momentum component'; ...
                           'units','meters seconds-1'; ...
                           'time','ocean_time' };                 
info.vars(end+1).name = 'temp';
info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                           'dimid',[0 4 8 11]; ...
                           'long_name','potential temperature'; ...
                           'units','Celsius'; ...
                           'time','ocean_time' };                      
info.vars(end+1).name = 'salt';
info.vars(end  ).info = {'type','NC_DOUBLE'; ...
                           'dimid',[0 4 8 11]; ...
                           'long_name','salinity'; ...
                           'units','Celsius'; ...
                           'time','ocean_time' };                        
                                               
% Last Info of initial condition variables; coordinates
if(spherical==1)
  info.vars(end-6).info(end+1,:) = {'coordinates','lon_rho, lat_rho, ocean_time'};
  info.vars(end-5).info(end+1,:) = {'coordinates','lon_u, lat_u, ocean_time'};
  info.vars(end-4).info(end+1,:) = {'coordinates','lon_v, lat_v, ocean_time'};
  info.vars(end-3).info(end+1,:) = {'coordinates','lon_u, lat_u, ocean_time'};
  info.vars(end-2).info(end+1,:) = {'coordinates','lon_v, lat_v, ocean_time'};
  info.vars(end-1).info(end+1,:) = {'coordinates','lon_rho, lat_rho, s_rho, ocean_time'};
  info.vars(end  ).info(end+1,:) = {'coordinates','lon_rho, lat_rho, s_rho, ocean_time'};
elseif(spherical==0)
  info.vars(end-6).info(end+1,:) = {'coordinates','x_rho, y_rho, ocean_time'};
  info.vars(end-5).info(end+1,:) = {'coordinates','x_u, y_u, ocean_time'};
  info.vars(end-4).info(end+1,:) = {'coordinates','x_v, y_v, ocean_time'};
  info.vars(end-3).info(end+1,:) = {'coordinates','x_u, y_u, ocean_time'};
  info.vars(end-2).info(end+1,:) = {'coordinates','x_v, y_v, ocean_time'};
  info.vars(end-1).info(end+1,:) = {'coordinates','x_rho, y_rho, s_rho, ocean_time'};
  info.vars(end  ).info(end+1,:) = {'coordinates','x_rho, y_rho, s_rho, ocean_time'};
end

% Done!
end