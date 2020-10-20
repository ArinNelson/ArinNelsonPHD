function [x,y,m] = r_grid_read(grid_file,grid_type,x_type)

  % Determine horizontal grid type
  switch x_type
    case 'x';	str_x = 'x';	str_y = 'y';
    case 'y';	str_x = 'lon';	str_y = 'lat';
    otherwise; error(['Unknown horizontal type: ' x_type]);
  end

  % Read in grid vars
  try x = ncread(grid_file,[x_type '_' grid_type]);    catch err; x=NaN; end
  try y = ncread(grid_file,[y_type '_' grid_type]);    catch err; y=NaN; end
  try m = ncread(grid_file,['mask_' grid_type]);       catch err; m=NaN; end

end