function [z,t] = roms_read_output_var(output_dir, output_prefix, var_name, depth_level)
%=========================================================================%
% [z, t] = roms_read_output_var(output_dir, output_prefix, var_name, depth_level)
% Gathers the variable 'var_name' at depth level 'depth_level' from the
% ROMS output residing in directory 'output_dir'.  The output files are
% assumed to have the format 'output_prefix'*.nc.
%=========================================================================%
% By Arin Nelson on 10/20/2020
% Last updated 02/08/2021
%=========================================================================%

% Ensure output_dir is in the proper format (ends in a \)
if( ~strcmp(output_dir(end),'\') ); output_dir = [output_dir '\'];  end

% Get the list of ROMS output files
file_list = ls([output_dir output_prefix '*.nc']);

% Get information on file and variable from first data file
file_info = ncinfo([output_dir file_list(1,:)]);
var_info  = file_info.Variables( strcmp({file_info.Variables.Name},var_name) );
if(isempty(var_info)); error(['Variable ' var_name ' not found in ROMS output!']); end

% Determine grid that variable resides on
switch var_info.Dimensions(1).Name(end)
  case 'o';   var_grid = 'rho';
  case 'u';   var_grid = 'u';
  case 'v';   var_grid = 'v';
  case 'i';   var_grid = 'psi';
end

% Determine if variable is 2D or 3D
var_dims = numel(var_info.Size)-1;

% Get time information
n_files     = size(file_list,1);
nt_per_file = var_info.Size(end);
nt          = n_files*nt_per_file;
  
% Initialize data variable
z = NaN(var_info.Size(1),var_info.Size(2),nt);

% If wanted, initialize time variable
if(nargout>1)
  t = NaN(nt,1);
end

% Read in data from files
n = 1;
for i=1:n_files
    
  % Indices
  ii = n:(n+nt_per_file-1);  
  
  % Load time
  if(nargout>1)
    t(ii) = ncread([output_dir file_list(i,:)],'ocean_time');
  end
  
  % Load data
  if(var_dims==2)
    tmp = ncread([output_dir file_list(i,:)],var_name);
  else
    tmp = squeeze( ncread([output_dir file_list(i,:)],var_name,[1 1 depth_level 1],[inf inf 1 inf]) );
  end
  
  % Save to array
  z(:,:,ii) = tmp;
  
  % Next step
  n = n + nt_per_file;

end
clear i n tmp;
