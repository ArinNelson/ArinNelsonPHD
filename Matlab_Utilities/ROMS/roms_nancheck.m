function results = roms_nancheck(ncfile)
% results = roms_nancheck(ncfile)
% checks the ROMS netcdf input file ncfile for any NaNs.  ROMS is in
% FORTRAN, which does not like NaNs.
%=========================================================================%
% By Arin Nelson on 03/10/2021
%=========================================================================%

  % Init output
  results = struct;

  % Get info from file
  info=ncinfo(ncfile);
  
  % Get list of variables
  varname = {info.Variables.Name};
  
  % Loop through variables and look for NaNs
  for i=1:numel(varname)
      
    % Load variable
    varvalue = ncread(ncfile,varname{i});
      
    % Look for NaNs
    ii = find(isnan(varvalue) | isinf(varvalue));
    if(~isempty(ii))
        results(end+1).varName = varname{i};
        results(end).nanIndex = ii;
    end
    clear varvalue ii;
      
  end
  clear i;

end