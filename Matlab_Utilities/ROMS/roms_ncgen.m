function status = roms_ncgen(fileName,fileInfo)
%=========================================================================%
% status = roms_ncgen(fileName,fileInfo)
% Create the NetCDF file 'fileName' from the information in the structure
% 'fileInfo'.  Typically fileInfo is generated by another file, such as
% roms_sample_output_var.m, etc., but you can use it for your own NetCDFs 
% given your fileInfo structure is constructed correctly.
%=========================================================================%
% By Arin Nelson on 09/23/2020
% Last updated 02/08/2021
%=========================================================================%

  % If all is successful, status will be 1
  status = 1;

  % Create the file
  ncid = netcdf.create(fileName,'CLOBBER');

  % Define the dimensions
  n_dim = size(fileInfo.dims,1);
  dimid = zeros(n_dim,1);
  for i=1:n_dim
    dimid(i) = netcdf.defDim(ncid,fileInfo.dims{i,1},fileInfo.dims{i,2});
  end
  clear i;

  % Put attributes
  n_att = size(fileInfo.atts,1);
  for i=1:n_att
    netcdf.putAtt(ncid,netcdf.getConstant('GLOBAL'),fileInfo.atts{i,1},fileInfo.atts{i,2});
  end

  % Define the variables
  for i=1:numel(fileInfo.vars)
    netcdf.defVar(ncid,fileInfo.vars(i).name,fileInfo.vars(i).info{1,2},fileInfo.vars(i).info{2,2});
  end

  % Put variable attributes
  for i=1:numel(fileInfo.vars)
    if(size(fileInfo.vars(i).info,1)>2)
      for j=3:size(fileInfo.vars(i).info,1)
        try
          netcdf.putAtt(ncid,i-1,fileInfo.vars(i).info{j,1},fileInfo.vars(i).info{j,2});
        catch err
          status = 0;
        end
      end
    end
  end

  % End definitions
  netcdf.endDef(ncid);

  % Close the file
  netcdf.close(ncid);

end