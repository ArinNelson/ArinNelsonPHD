function status = r_ncgen_std_ini(ini_file, grd_file, out_file)
%=========================================================================%
% status = r_ncgen_std_ini(infile, grd_file, outfile)
% Generates a template error covariance of initial conditions file 
% <out_file> using values as found in the initial conditions data file
% <ini_file> and grid file <grd_file>.
%=========================================================================%
% by Arin Nelson on 09/26/2020
%=========================================================================%

% Information on dimensions, etc. can be gathered using ncinfo
ini_info = ncinfo(ini_file);
grd_info = ncinfo(grd_file);

% The new variable dimensions are essentially the same as the ini file's
out_info = struct;
out_info.Dimensions = ini_info.Dimensions;

% But the global attributes are set separately
out_info.Attributes = [];
out_info.Attributes(1).Name  = 'type';
out_info.Attributes(1).Value = '4D-Var error covariance standard deviation';
out_info.Attributes(2).Name  = 'Conventions';
out_info.Attributes(2).Value = 'CF-1.4';
out_info.Attributes(3).Name = 'grd_file';
out_info.Attributes(3).Value = grd_file;
out_info.Attributes(4).Name  = 'history';
out_info.Attributes(4).Value = ['File generated using r_ncgen_std_ini on ' datestr(now)];

% All of the variables from the ini file are used
out_info.Variables = ini_info.Variables;

% And some variables from the grid file are used as well
grd_vars = {'angle','mask_rho','mask_u','mask_v'};
for i=1:numel(grd_vars)
  out_info.Variables(end+1) = grd_info.Variables(  strcmp({grd_info.Variables.Name},grd_vars{i})==1  );
end

% Some other things that have to be set
out_info.Name = '/';
out_info.Format = ini_info.Format;

% Generate the std file
ncwriteschema(out_file,out_info);

% Write variables present in ini file to this file
ini_vars = {'spherical','Vtransform','Vstretching','theta_s','theta_b','Tcline','hc','s_rho','s_w','Cs_r','Cs_w',...
            'lon_rho','lat_rho','lon_u','lat_u','lon_v','lat_v','ocean_time'};
for i=1:numel(ini_vars)
  ncwrite(out_file,ini_vars{i},ncread(ini_file,ini_vars{i}));
end

% Write variables present in grid file to this file
grd_vars = {'angle','mask_rho','mask_u','mask_v'};
for i=1:numel(grd_vars)
  ncwrite(out_file,grd_vars{i},ncread(grd_file,grd_vars{i}));
end

% Lastly, write some default values to the momentum and tracer variables
% But first, we need the land/sea masks
mask_rho = ncread(grd_file,'mask_rho');
mask_u   = ncread(grd_file,'mask_u'  );
mask_v   = ncread(grd_file,'mask_v'  );

% Now, generate and write the default values
N  = out_info.Dimensions(  strcmp({out_info.Dimensions.Name},'s_rho')==1  ).Length;
zeta = zeros(size(mask_rho));	
  zeta(mask_rho==0) = 0; 
  zeta(mask_rho==1) = 0.01;    
  zeta = repmat(zeta,[1 1 1]);        
  ncwrite(out_file,'zeta',zeta);   clear zeta;
ubar = zeros(size(mask_u));     
  ubar(mask_u==0) = 0;    
  ubar(mask_u==1) = 0.01;    
  ubar = repmat(ubar,[1 1 1]);        
  ncwrite(out_file,'ubar',ubar);   clear ubar;
vbar = zeros(size(mask_v));     
  vbar(mask_v==0) = 0;  
  vbar(mask_v==1) = 0.01;    
  vbar = repmat(vbar,[1 1 1]);        
  ncwrite(out_file,'vbar',vbar);   clear vbar;
u = zeros(size(mask_u));     
  u(mask_u==0) = 0;    
  u(mask_u==1) = 0.01;       
  u = repmat(u,[1 1 N 1]);     
  ncwrite(out_file,'u',u);         clear u;
v = zeros(size(mask_v));
  v(mask_v==0) = 0;    
  v(mask_v==1) = 0.01;       
  v = repmat(v,[1 1 N 1]);     
  ncwrite(out_file,'v',v);         clear v;
temp = zeros(size(mask_rho));     
  temp(mask_rho==0) = 0;	
  temp(mask_rho==1) = 0.01;	
  temp = repmat(temp,[1 1 N 1]);     
  ncwrite(out_file,'temp',temp);   clear temp;
salt = zeros(size(mask_rho));     
  salt(mask_rho==0) = 0;	
  salt(mask_rho==1) = 0.01;	
  salt = repmat(salt,[1 1 N 1]);     
  ncwrite(out_file,'salt',salt);   clear temp;  
clear mask_*  

% If made it this far, code was successful!
status = 1;

end %function