function [t, nt_per_file] = roms_read_output_times(output_dir, output_prefix)
%=========================================================================%
% t = roms_read_output_times(output_dir, output_prefix)
% Gathers the variable 'ocean_time' from the ROMS output residing in 
% directory 'output_dir'.  The output files are assumed to have the format 
% 'output_prefix'*.nc.
%=========================================================================%
% By Arin Nelson on 02/08/2021
%=========================================================================%

  % Ensure output_dir is in the proper format (ends in a \)
  if( ~strcmp(output_dir(end),'\') ); output_dir = [output_dir '\'];  end

  % Get the list of ROMS output files
  file_list = ls([output_dir output_prefix '*.nc']);

  % Get number of time steps per file
  nt_per_file = numel(ncread([output_dir file_list(1,:)],'ocean_time'));

  % Get time information
  n_files     = size(file_list,1);
  nt          = n_files*nt_per_file;
  t           = NaN(nt,1);

  % Read in ocean_time from each file
  n = 1;
  for i=1:n_files
    ii = n:(n+nt_per_file-1);  
    t(ii) = ncread([output_dir file_list(i,:)],'ocean_time');
    n = n + nt_per_file;
  end
  
  % Clean-up
  clear file_list nt_per_file n_files nt n ii i;

end