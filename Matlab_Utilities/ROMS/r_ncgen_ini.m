function status = r_ncgen_ini(in_file, out_file)
%=========================================================================%
% status = r_ncgen_ini(infile, outfile)
% Generates a template initial conditions file <out_file> using values as 
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
out_info.dims(  strcmp({out_info.dims.name},'xi_psi' ) ==1  ).value = Lm-1;
out_info.dims(  strcmp({out_info.dims.name},'eta_rho') ==1  ).value = Mm;
out_info.dims(  strcmp({out_info.dims.name},'eta_u'  ) ==1  ).value = Mm;
out_info.dims(  strcmp({out_info.dims.name},'eta_v'  ) ==1  ).value = Mm-1;
out_info.dims(  strcmp({out_info.dims.name},'eta_psi') ==1  ).value = Mm-1;
out_info.dims(  strcmp({out_info.dims.name},'s_rho'  ) ==1  ).value = N;
out_info.dims(  strcmp({out_info.dims.name},'s_w'    ) ==1  ).value = N+1;

% Other dimension lengths (# active tracers, time steps)
% NOTE: May need to change tracer dimension calculation to include passive, sediment, and other tracers!
out_info.dims(  strcmp({out_info.dims.name},'tracer'    )==1  ).value = str2double(  in_vars(  strcmp({in_vars.name},'NAT')==1  ).value{1}  );
out_info.dims(  strcmp({out_info.dims.name},'ocean_time')==1  ).value = netcdf.getConstant('UNLIMITED');

% With the dimensions set, generate the initial conditions file
r_ncgen(out_file,out_info);

% The horizontal grid variables can be transferred from the grid file
grid_vars = {'spherical','lon_rho','lat_rho','lon_u','lat_u','lon_v','lat_v'};
for i=1:numel(grid_vars)
  ncwrite(out_file,grid_vars{i},ncread(grdname,grid_vars{i}));
end
clear i grid_vars;

% Some vertical grid variables can be read in from the roms.in file
Vtransform  = str2double(  in_vars(  strcmp({in_vars.name},'Vtransform' )==1  ).value{1}  );
Vstretching = str2double(  in_vars(  strcmp({in_vars.name},'Vstretching')==1  ).value{1}  );
theta_s     = str2double(  in_vars(  strcmp({in_vars.name},'THETA_S'    )==1  ).value{1}  );
theta_b     = str2double(  in_vars(  strcmp({in_vars.name},'THETA_B'    )==1  ).value{1}  );
Tcline      = str2double(  in_vars(  strcmp({in_vars.name},'TCLINE'     )==1  ).value{1}  );
hc          = Tcline;

% Others have to be computed depending on these previous values
[s_rho,Cs_r] = stretching(Vstretching,theta_s,theta_b,hc,N,0,0);
[s_w,  Cs_w] = stretching(Vstretching,theta_s,theta_b,hc,N,1,0);

% Now the vertical grid variables can be written to the i.c. NetCDF file
vert_vars = {'Vtransform','Vstretching','theta_s','theta_b','Tcline','hc','s_rho','s_w','Cs_r','Cs_w'};
for i=1:numel(vert_vars)
  eval(['ncwrite(out_file,vert_vars{i},' vert_vars{i} ');']);  
end
clear i vert_vars

% Lastly, write some default values to the momentum and tracer variables
% But first, we need the land/sea masks
mask_rho = ncread(grdname,'mask_rho');
mask_u   = ncread(grdname,'mask_u'  );
mask_v   = ncread(grdname,'mask_v'  );

% Now, generate and write the default values
zeta = zeros(size(mask_rho));	
  zeta(mask_rho==0) = realmax('single');    
  zeta(mask_rho==1) = 0;    
  zeta = repmat(zeta,[1 1 1]);        
  ncwrite(out_file,'zeta',zeta);   clear zeta;
ubar = zeros(size(mask_u));     
  ubar(mask_u==0) = realmax('single');    
  ubar(mask_u==1) = 0;    
  ubar = repmat(ubar,[1 1 1]);        
  ncwrite(out_file,'ubar',ubar);   clear ubar;
vbar = zeros(size(mask_v));     
  vbar(mask_v==0) = realmax('single');    
  vbar(mask_v==1) = 0;    
  vbar = repmat(vbar,[1 1 1]);        
  ncwrite(out_file,'vbar',vbar);   clear vbar;
u = zeros(size(mask_u));     
  u(mask_u==0) = realmax('single');    
  u(mask_u==1) = 0;       
  u = repmat(u,[1 1 N 1]);     
  ncwrite(out_file,'u',u);         clear u;
v = zeros(size(mask_v));
  v(mask_v==0) = realmax('single');    
  v(mask_v==1) = 0;       
  v = repmat(v,[1 1 N 1]);     
  ncwrite(out_file,'v',v);         clear v;
temp = zeros(size(mask_rho));     
  temp(mask_rho==0) = realmax('single');	
  temp(mask_rho==1) = 20;	
  temp = repmat(temp,[1 1 N 1]);     
  ncwrite(out_file,'temp',temp);   clear temp;
salt = zeros(size(mask_rho));     
  salt(mask_rho==0) = realmax('single');	
  salt(mask_rho==1) = 20;	
  salt = repmat(salt,[1 1 N 1]);     
  ncwrite(out_file,'salt',salt);   clear temp;  
clear mask_*  

% If made it this far, code was successful!
status = 1;

end %function