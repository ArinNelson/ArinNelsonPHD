function status = r_ncgen_bry(in_file, out_file)
%=========================================================================%
% status = r_ncgen_bry(infile, outfile)
% Generates a template boundary conditions file <out_file> using values as 
% found in the roms.in file <in_file>
%=========================================================================%
% by Arin Nelson on 09/26/2020
%=========================================================================%

% Check if roms.in file exists.  If it does, read it in to memory
if(exist(in_file,'file')~=2)
  error(['Cannot find specified roms.in file: ' in_file]);
else
  in_vars = r_in_read(in_file);
end

% Get spherical switch from grid file
grdname   = in_vars(  strcmp({in_vars.name},'GRDNAME')==1  ).value{1};
spherical = ncread(grdname,'spherical');

% Generate template of initial conditions NetCDF file info 
out_info = r_nctemplate_ini(spherical);

% Get grid dimensions from roms.in
Lm = str2double(  in_vars(  strcmp({in_vars.name},'Lm')==1  ).value{1}  ) + 2; 
Mm = str2double(  in_vars(  strcmp({in_vars.name},'Mm')==1  ).value{1}  ) + 2;
N  = str2double(  in_vars(  strcmp({in_vars.name},'N' )==1  ).value{1}  );

% Set grid-related dimension lengths to initial conditions NetCDF file info
out_info.dims(  strcmp({out_info.dims.name},'xi_rho' ) ==1  ).value = Lm;
out_info.dims(  strcmp({out_info.dims.name},'xi_u'   ) ==1  ).value = Lm-1;
out_info.dims(  strcmp({out_info.dims.name},'xi_v'   ) ==1  ).value = Lm;
out_info.dims(  strcmp({out_info.dims.name},'eta_rho') ==1  ).value = Mm;
out_info.dims(  strcmp({out_info.dims.name},'eta_u'  ) ==1  ).value = Mm;
out_info.dims(  strcmp({out_info.dims.name},'eta_v'  ) ==1  ).value = Mm-1;
out_info.dims(  strcmp({out_info.dims.name},'s_rho'  ) ==1  ).value = N;
out_info.dims(  strcmp({out_info.dims.name},'s_w'    ) ==1  ).value = N+1;

% Other dimension lengths (# active tracers, time steps)
% NOTE: May need to change tracer dimension calculation to include passive, sediment, and other tracers!
out_info.dims(  strcmp({out_info.dims.name},'tracer'  )==1  ).value = str2double(  in_vars(  strcmp({in_vars.name},'NAT')==1  ).value{1}  );
out_info.dims(  strcmp({out_info.dims.name},'bry_time')==1  ).value = netcdf.getConstant('UNLIMITED');

% With the dimensions set, generate the initial conditions file
r_ncgen(out_file,out_info);

% HERE









end %function