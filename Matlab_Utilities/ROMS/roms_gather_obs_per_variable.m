function status = roms_gather_obs_per_variable(obs_info,mdl_info,options)
% status = roms_gather_obs_per_variable(obs_info,mdl_info,options
% Reads observations from the information given by structure 'obs_info',
% translates things based on model information in 'mdl_info', and generates
% a single file per variable.
% 
% Required pieces of obs_info are:
% 
%  FileName (1xN cell array of strings)
%     cell array of strings containing path to observation files
% 
%  VarName (1xN cell array of cells)
%     cell array of cells of strings listing the variables in the 
%     observation files in order of the ROMS variable list:
%     [1 2 3 4 5 6 7 (...)]={zeta, ubar, vbar, u, v, temp, salt, (others)}
%     Note: if any are unavailable per file, keep that cell entry empty.
% 
%  SampleTime (1xN cell array of 1D vectors)
%     cell array of vectors containing the sampling times relative to the
%     observations in each file in units of days from the reference time in
%     the model's roms.in file.
% 
% 
% One of the following:
% 
%  Lon or XI (1xN numerical array)
%    XI is the fractional x-grid location of the observation.  If Lon is
%    provided instead, XI is computed from the grid file specified in
%    'mdl_info'.
%    Note: xi's of rho-grid cell centers are 0.5 : 1 : (L-0.5)
% 
%  Lat or YI (1xN numerical array)
%    YI is the fractional y-grid location of the observation.  If Lat is
%    provided instead, YI is computed from the grid file specified in
%    'mdl_info'.
%    Note: yi's of rho-grid cell centers are 0.5 : 1 : (M-0.5)
% 
%  Depth or ZI (1xN cell array)
%    ZI is the fractional z-level location(s) of the observations.
%    Depth can be given instead (ROMS will translate to ZI internally).
%    Note: ROMS lists depths in reverse order, so ZI=1 is the bottom level,
%          and ZI=(max) is the surface level.
% 
% 
% Required pieces in mdl_info are:
%  
%  GridFile (string)
%     directory to the grid file used by the model
%
%  
% 
%==========================================================================
% by Arin Nelson
% on 02/10/2021
%==========================================================================
