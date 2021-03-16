function varargout = roms_read_outputs_at_time(output_dir, output_prefix, var_names, depth_level, sample_index)
%=========================================================================%
% [z, t] = roms_read_output_var(output_dir, output_prefix, var_names, depth_level, time)
% Gathers the variables 'var_names' at depth level 'depth_level' from the
% ROMS output residing in directory 'output_dir' at time 'sample_time'.  
% The output files are assumed to have the format 'output_prefix'*.nc.
%=========================================================================%
% By Arin Nelson on 02/08/2021
%=========================================================================%

% Ensure output_dir is in the proper format (ends in a \)
if( ~strcmp(output_dir(end),'\') ); output_dir = [output_dir '\'];  end

% Get the list of ROMS output files
file_list = ls([output_dir output_prefix '*.nc']);

% Number of variables
n_var = numel(var_names);

% Get information on file and variable from first data file
file_info = ncinfo([output_dir file_list(1,:)]);
var_info  = cell(n_var,1);
for i=1:n_var
  var_info{i} = file_info.Variables( strcmp({file_info.Variables.Name},var_names{i}) );
  if(isempty(var_info{i})); error(['Variable ' var_names{i} ' not found in ROMS output!']); end
end

% Determine grid that variable resides on
var_grid = cell(n_var,1);
for i=1:n_var
  switch var_info{i}.Dimensions(1).Name(end)
    case 'o';   var_grid{i} = 'rho';
    case 'u';   var_grid{i} = 'u';
    case 'v';   var_grid{i} = 'v';
    case 'i';   var_grid{i} = 'psi';
  end
end

% Determine if variable is 2D or 3D
var_dims = zeros(n_var,1);
for i=1:n_var
  var_dims(i) = numel(var_info{i}.Size)-1;
end

% Get time information
n_files     = size(file_list,1);
nt_per_file = var_info{1}.Size(end);
nt          = n_files*nt_per_file;

% Find file associated with wanted time
ifile = ceil(sample_index/nt_per_file);
isamp = mod(sample_index,nt_per_file); if(isamp==0); isamp=nt_per_file; end

% Read in variables at this time
varargout = cell(n_var,1);
for i=1:n_var

  % Load data
  if(var_dims==2)
    varargout{i} = ncread([output_dir file_list(ifile,:)],var_names{i},[1 1 isamp],[inf inf 1]);
  else
    varargout{i} = squeeze( ncread([output_dir file_list(i,:)],var_names{i},[1 1 depth_level isamp],[inf inf 1 1]) );
  end

end
clear i;

end