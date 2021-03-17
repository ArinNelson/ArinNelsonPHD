% Example: setting up an operation to take synthetic observations from a
% ROMS output.  Another code will translate the observations into a DA
% observations input file.
% by Arin Nelson
% on 02/08/2021
%=========================================================================%

%--------------------------------------------------------------------------
% Example #1: multiple platforms, surface and bottom only, same variables
% each

% Observations Info
% N. Prudence, Conimicut Pt., Mtn. View, Quonset Point, Mt. Hope Bay, 
% Poppasquash Pt., and Sally Rock

% obsInfo.Name  = {'NP','CP','MV','QP','MHB','PP','SR'};
% obsInfo.Label = {'NP','CP','MV','QP','MHB','PP','SR'};
% obsInfo.Lon   = -71 - [3549 3437 3909 3800 2156 3181 4243]./10000;
% obsInfo.Lat   =  41 + [6708 7128 6385 5903 6799 6492 6760]./10000;
% obsInfo.Time  = 0:(1/12):(60-(1/12));
% obsInfo.TI    = 1:1:720;
% obsInfo.Vars  = 'zeta,temp,salt';
% obsInfo.Depth = [15 1];
% obsInfo.ZI    = [15 1];
% 
% % Models Info
% mdlInfo.Name  = 'NGBay_Case171';             % Simulation name
% mdlInfo.Label = 'ngbay_case171';
% mdlInfo.Dir   = 'F:\Datasets\ROMS\case171\'; % Directory containing the simulation's outputs
% mdlInfo.Prfx  = 'ocean_avg_0';              % All output files have this format
% mdlInfo.Grid  = 'C:\Library\ROMS_Stuff\Resources\ngbay_grd.nc';  % Grid file associated with this simulation
% 
% % Compute XI and YI
% [mdlx,mdly,mdlm] = roms_read_grid(mdlInfo.Grid,'rho','lon');
% mdli=1:1000;    %0.5:(1000-0.5);
% mdlj=1:1100;    %0.5:(1100-0.5);
% [ii,jj] = meshgrid(mdlj,mdli);
% kk = ii + sqrt(-1).*jj;
% ntrplnt = scatteredInterpolant(mdlx(mdlm==1),mdly(mdlm==1),kk(mdlm==1),'linear');
% lm = ntrplnt(obsInfo.Lon,obsInfo.Lat);
% obsInfo.XI   = real(lm);
% obsInfo.YI   = imag(lm);
% %test = [mdlx(round(obsInfo.XI),round(obsInfo.YI)), mdly(round(obsInfo.XI),round(obsInfo.YI))] - [obsInfo.Lon, obsInfo.Lat];

% % Perform
% roms_sample_output(mdlInfo,obsInfo);

%--------------------------------------------------------------------------
% Example #2: combining observation files into a single DA-OBS file

% Observation info
nObs             = 7;
obsInfo.obsLabel = {'NP','CP','MV','QP','MHB','PP','SR'};
obsInfo.obsFile  = cell(nObs,1);    for i=1:numel(obsInfo.obsLabel);    obsInfo.obsFile{i} = ['ngbay_case171_obs_' obsInfo.obsLabel{i} '.nc'];  end
obsInfo.obsXI    = zeros(nObs,1);
obsInfo.obsYI    = zeros(nObs,1);
obsInfo.obsDepth = [15 1];
obsInfo.obsZI    = [15 1];
obsInfo.obsTime  = ncread(obsInfo.obsFile{1},'time');
obsInfo.obsVars  = {'zeta','temp','salt'};

% Model info
mdlInfo.Title = 'NGBay_Case171';             % Simulation name
mdlInfo.Label = 'ngbay_case171';
mdlInfo.Dir   = 'F:\Datasets\ROMS\case171\'; % Directory containing the simulation's outputs
mdlInfo.Prfx  = 'ocean_avg_0';              % All output files have this format
mdlInfo.Grid  = 'C:\Library\ROMS_Stuff\Resources\ngbay_grd.nc';  % Grid file associated with this simulation
mdlInfo.LM    = [330 360];

% Do translation
roms_translate_obs_to_da(obsInfo,mdlInfo)