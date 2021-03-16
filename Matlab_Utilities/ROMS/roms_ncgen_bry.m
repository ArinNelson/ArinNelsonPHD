function status = roms_ncgen_bry(in_file,varargin)
%=========================================================================%
% status = roms_ncgen_bry(infile)
% Generates a template boundary conditions file out_file using values as 
% found in the roms.in file in_file
%=========================================================================%
% by Arin Nelson on 09/26/2020
% last updated 02/27/2021
%=========================================================================%

% NOTE:
% W: (1,:)
% S: (:,1)
% E: (end,:)
% N: (:,end)

% Initial stuff
strWSEN = {'west','south','east','north'};

% Check if roms.in file exists.  If it does, read it in to memory
if(exist(in_file,'file')~=2)
  error(['Cannot find specified roms.in file: ' in_file]);
else
  in_vars = roms_read_infile(in_file);
end

% Get spherical switch from grid file
grid_file   = in_vars(  strcmp({in_vars.name},'GRDNAME')==1  ).value{1};
spherical = ncread(grid_file,'spherical');

% Get WSEN boundary info from grid file
bryinfo = in_vars( strcmp({in_vars.name},'LBC(isFsur)')==1 ).value(1:4);
WSEN = ones(4,1);
for i=1:4
  if( strcmp(bryinfo{i},'Clo') ); WSEN(i)=0; end
end

% Generate template of boundary conditions NetCDF file info 
out_info = roms_nctemplate_bry(spherical,WSEN);

% Get grid dimensions from roms.in
Lm = in_vars(  strcmp({in_vars.name},'Lm')==1  ).value{1} + 2; 
Mm = in_vars(  strcmp({in_vars.name},'Mm')==1  ).value{1} + 2;
N  = in_vars(  strcmp({in_vars.name},'N' )==1  ).value{1};

% Set grid-related dimension lengths to initial conditions NetCDF file info
out_info.dims{ strcmp(out_info.dims(:,1),'xi_rho' ) ==1  , 2} = Lm;
out_info.dims{ strcmp(out_info.dims(:,1),'xi_u'   ) ==1  , 2} = Lm-1;
out_info.dims{ strcmp(out_info.dims(:,1),'xi_v'   ) ==1  , 2} = Lm;
%out_info.dims{ strcmp(out_info.dims(:,1),'xi_psi' ) ==1  , 2} = Lm-1;
out_info.dims{ strcmp(out_info.dims(:,1),'eta_rho') ==1  , 2} = Mm;
out_info.dims{ strcmp(out_info.dims(:,1),'eta_u'  ) ==1  , 2} = Mm;
out_info.dims{ strcmp(out_info.dims(:,1),'eta_v'  ) ==1  , 2} = Mm-1;
%out_info.dims{ strcmp(out_info.dims(:,1),'eta_psi') ==1  , 2} = Mm-1;
out_info.dims{ strcmp(out_info.dims(:,1),'s_rho'  ) ==1  , 2} = N;
out_info.dims{ strcmp(out_info.dims(:,1),'s_w'    ) ==1  , 2} = N+1;

% Other dimension lengths (# active tracers, time steps)
% NOTE: May need to change tracer dimension calculation to include passive, sediment, and other tracers!
out_info.dims{ strcmp(out_info.dims(:,1),'tracer'  )==1 , 2} = in_vars(  strcmp({in_vars.name},'NAT')==1  ).value{1};
out_info.dims{ strcmp(out_info.dims(:,1),'bry_time')==1 , 2} = netcdf.getConstant('UNLIMITED');

% Get bry file name from roms.in
out_file = in_vars( strcmp({in_vars.name},'BRYNAME')==1 ).value{1};

% With the dimensions set, generate the initial conditions file
roms_ncgen(out_file,out_info);

% Grid variables to write out depend on WSEN
grid_vars = {'lon_rho','lat_rho','lon_u','lat_u','lon_v','lat_v'};
if(WSEN(1)==1)
  for j=1:numel(grid_vars)
    tmp = ncread(grid_file,grid_vars{j});
    ncwrite(out_file,[grid_vars{j} '_west'],tmp(1,:));
  end
  clear j tmp;
end
if(WSEN(2)==1)
  for j=1:numel(grid_vars)
    tmp = ncread(grid_file,grid_vars{j});
    ncwrite(out_file,[grid_vars{j} '_south'],tmp(:,1));
  end
  clear j tmp;
end
if(WSEN(3)==1)
  for j=1:numel(grid_vars)
    tmp = ncread(grid_file,grid_vars{j});
    ncwrite(out_file,[grid_vars{j} '_east'],tmp(end,:));
  end
  clear j tmp;
end
if(WSEN(4)==1)
  for j=1:numel(grid_vars)
    tmp = ncread(grid_file,grid_vars{j});
    ncwrite(out_file,[grid_vars{j} '_north'],tmp(:,end));
  end
  clear j tmp;
end

% Some vertical grid variables can be read in from the roms.in file
Vtransform  = in_vars( strcmp({in_vars.name},'Vtransform' )==1  ).value{1};
Vstretching = in_vars( strcmp({in_vars.name},'Vstretching')==1  ).value{1};
theta_s     = in_vars( strcmp({in_vars.name},'THETA_S'    )==1  ).value{1};
theta_b     = in_vars( strcmp({in_vars.name},'THETA_B'    )==1  ).value{1};
Tcline      = in_vars( strcmp({in_vars.name},'TCLINE'     )==1  ).value{1};
hc          = Tcline;

% Others have to be computed depending on these previous values
[s_rho,Cs_r] = stretching(Vstretching,theta_s,theta_b,hc,N,0,0);
[s_w,  Cs_w] = stretching(Vstretching,theta_s,theta_b,hc,N,1,0);

% Now the vertical grid variables can be written to the b.c. NetCDF file
vert_vars = {'Vtransform','Vstretching','theta_s','theta_b','Tcline','hc','s_rho','s_w','Cs_r','Cs_w'};
for i=1:numel(vert_vars)
  try eval(['ncwrite(out_file,vert_vars{i},' vert_vars{i} ');']);  
  catch err
      pause(1e-9);
  end
end
clear i vert_vars

% Lastly, generate some default values?
% NOTE: This is not straightforward as it depends on the time range the
% user will run the simulation for, and a user-specified dt
% RIROMS bry forcing currently uses dt=6hrs, with t=0:21600:31536000 (units
% of seconds from reference time in roms.in)
if(~isempty(varargin))
  bry_time = varargin{1};  
  nt = numel(bry_time);  
  ncwrite(out_file,'bry_time',bry_time);

  % West & East variables
  if(WSEN(1)==1)
    ncwrite(out_file,'zeta_west',zeros(Mm,nt));
    ncwrite(out_file,'ubar_west',zeros(Mm,nt));
    ncwrite(out_file,'vbar_west',zeros(Mm-1,nt));
    ncwrite(out_file,'u_west',zeros(Mm,N,nt));
    ncwrite(out_file,'v_west',zeros(Mm-1,N,nt));
    ncwrite(out_file,'temp_west',ones(Mm,N,nt).*20);
    ncwrite(out_file,'salt_west',ones(Mm,N,nt).*30);
  end
  if(WSEN(3)==1)
    ncwrite(out_file,'zeta_east',zeros(Mm,nt));
    ncwrite(out_file,'ubar_east',zeros(Mm,nt));
    ncwrite(out_file,'vbar_east',zeros(Mm-1,nt));
    ncwrite(out_file,'u_east',zeros(Mm,N,nt));
    ncwrite(out_file,'v_east',zeros(Mm-1,N,nt));
    ncwrite(out_file,'temp_east',ones(Mm,N,nt).*20);
    ncwrite(out_file,'salt_east',ones(Mm,N,nt).*30);
  end
  
  % South & North variables
  if(WSEN(2)==1)
    ncwrite(out_file,'zeta_south',zeros(Lm,nt));
    ncwrite(out_file,'ubar_south',zeros(Lm-1,nt));
    ncwrite(out_file,'vbar_south',zeros(Lm,nt));
    ncwrite(out_file,'u_south',zeros(Lm-1,N,nt));
    ncwrite(out_file,'v_south',zeros(Lm,N,nt));
    ncwrite(out_file,'temp_south',ones(Lm,N,nt).*20);
    ncwrite(out_file,'salt_south',ones(Lm,N,nt).*30);
  end  
  if(WSEN(4)==1)
    ncwrite(out_file,'zeta_north',zeros(Lm,nt));
    ncwrite(out_file,'ubar_north',zeros(Lm-1,nt));
    ncwrite(out_file,'vbar_north',zeros(Lm,nt));
    ncwrite(out_file,'u_north',zeros(Lm-1,N,nt));
    ncwrite(out_file,'v_north',zeros(Lm,N,nt));
    ncwrite(out_file,'temp_north',ones(Lm,N,nt).*20);
    ncwrite(out_file,'salt_north',ones(Lm,N,nt).*30);
  end 
  
end

end %function