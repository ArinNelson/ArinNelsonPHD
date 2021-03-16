function status = roms_ncgen_frc(in_file,var_list,out_file)
%=========================================================================%
% status = roms_ncgen_ini(infile,varlist,outfile)
% Generates a template surface forcings file outfile using values as found 
% in the roms.in file in_file and the variables present in cell array
% var_list.
%=========================================================================%
% by Arin Nelson on 03/08/2021
%=========================================================================%

  % Check if roms.in file exists.  If it does, read it in to memory
  if(exist(in_file,'file')~=2)
    error(['Cannot find specified roms.in file: ' in_file]);
  else
    in_vars = roms_read_infile(in_file);
  end

  % Get spherical switch from grid file
  grdname   = in_vars(  strcmp({in_vars.name},'GRDNAME')==1  ).value{1};
  spherical = ncread(grdname,'spherical');

  % Generate template of initial conditions NetCDF file info 
  if(strcmp(spherical,'T')==1 || spherical==1); spherical=1; else; spherical=0; end
  out_info = roms_nctemplate_frc(var_list);

  % Get grid dimensions from roms.in
  %Lm = in_vars(  strcmp({in_vars.name},'Lm')==1  ).value{1} + 2; 
  %Mm = in_vars(  strcmp({in_vars.name},'Mm')==1  ).value{1} + 2;
  Lm = 2;
  Mm = 2;

  % Set grid-related dimension lengths to initial conditions NetCDF file info
  out_info.dims{  strcmp( out_info.dims(:,1),'nlon'  ) ==1  , 2} = Lm;
  out_info.dims{  strcmp( out_info.dims(:,1),'nlat' ) ==1  , 2} = Mm;
  out_info.dims{  strcmp( out_info.dims(:,1),'time' ) ==1, 2} = 0;
  
  % Set time-related attributes
  for i=1:numel(var_list)
    out_info.vars(i+2).info{4,2} = 'days since 2006/01/01 00:00:00';  
  end
  
  % With the dimensions set, generate the initial conditions file
  roms_ncgen(out_file,out_info);
  
  % The horizontal grid variables can be transferred from the grid file
  tmpx = ncread(grdname,'lon_rho');     tmpx = tmpx([1 end],[1 end]) + [[-1 -1];[1 1]]./100;   ncwrite(out_file,'lon',tmpx);
  tmpy = ncread(grdname,'lat_rho');     tmpy = tmpy([1 end],[1 end]) + [[-1 1];[-1 1]]./100;   ncwrite(out_file,'lat',tmpy);
  
  %ncwrite(out_file,'lon',ncread(grdname,'lon_rho'));
  %ncwrite(out_file,'lat',ncread(grdname,'lat_rho'));
  
  % THAT'S IT!  NO DEFAULT VALUES THOUGH...