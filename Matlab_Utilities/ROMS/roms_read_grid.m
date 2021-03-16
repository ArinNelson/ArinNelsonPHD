function [x, y, m] = roms_read_grid(grid_file, grid_type, x_type)
%=========================================================================%
% [x,y,m] = roms_read_grid(grid_file, grid_type, x_type)
% Reads in x, y, and mask data from the ROMS grid NetCDF file 'grid_file'.
% 'grid_type' can be 'rho', 'u', 'v', or 'psi'.
% 'x_type' can be 'x' or 'lon'.
%=========================================================================%
% By Arin Nelson on 09/23/2020
% Last updated 02/08/2021
%=========================================================================%

% SOME ADDITIONAL NOTES
% GRID CELL AREA IS 1/(pm.*pm) (RHO GRID ONLY)

  % Determine horizontal grid type
  switch x_type
    case 'x';	str_x = 'x';	str_y = 'y';
    case 'lon';	str_x = 'lon';	str_y = 'lat';
    otherwise; error(['Unknown horizontal type: ' x_type]);
  end

  % Read in grid vars
  try x = ncread(grid_file,[str_x '_' grid_type]);    catch err; x=NaN; end
  try y = ncread(grid_file,[str_y '_' grid_type]);    catch err; y=NaN; end
  try m = ncread(grid_file,['mask_' grid_type]);      catch err; m=NaN; end

end