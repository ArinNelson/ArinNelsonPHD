% Example/How-To: Checking River Positions and Outflow Directions
% Arin Nelson
% 03/16/2021
%=========================================================================%
clc; close all;

% % Needed files
% grdFile = 'ngbay_grd.nc';
% rvrFile = 'C:\Library\ROMS_Stuff\Resources\forcefiles_riroms\riroms_rivers_2006-2010_dt24hrs.nc';

% % Needed variables
% mr = ncread(grdFile,'mask_rho');
% xr = 0.5:size(mr,1);
% yr = 0.5:size(mr,2);
% rx = ncread(rvrFile,'river_Xposition');
% ry = ncread(rvrFile,'river_Eposition');
% rd = ncread(rvrFile,'river_direction');
% rt = nanmean(ncread(rvrFile,'river_transport'),2);
% nr = numel(rx);

% % Show u rivers
% i = find(rd==0);
% figure; hold on;
%   imagesc(xr-0.5,yr-0.5,mr'); set(gca,'ydir','normal'); axis tight; box on;
%   plot(rx(i)-0.5,ry(i),'or');
%   quiver(rx(i)-0.5,ry(i),sign(rt(i)),zeros(numel(i),1),1e-2,'-k');
%   for j=i; text(rx(j)-0.5,ry(j)+0.25,num2str(j)); end
% hold off;

% % Show v rivers
% i = find(rd==1);
% figure; hold on;
%   imagesc(xr-0.5,yr-0.5,mr'); set(gca,'ydir','normal'); axis tight; box on;
%   plot(rx(i),ry(i)-0.5,'or');
%   quiver(rx(i),ry(i)-0.5,zeros(numel(i),1),sign(rt(i)),1e-2,'-k');
%   for j=i; text(rx(j),ry(j)-0.25,num2str(j)); end
% hold off;

% My new grid
newFile = 'C:\Library\ROMS_Stuff\Resources\riroms_grd_thirded_v2_trunc.nc';
rx_new  = round(rx.*(330/1000));
ry_new  = round(ry.*(360/1100));
rt_new  = rt;
rd_new  = rd;
mr_new  = ncread(newFile,'mask_rho');
xr_new  = 0.5:size(mr_new,1);
yr_new  = 0.5:size(mr_new,2);

% Fixes
rx_new(03) = 193;               ry_new(03) = 300;       rd_new(03) = 0;     rt_new(03) = 1-rt(03);  % NOTE SIGN CHANGE!
rx_new(04) = 193;               ry_new(04) = 301;       rd_new(04) = 0;     rt_new(04) = 1-rt(04);  % NOTE SIGN CHANGE!
rx_new(05) = 190;               ry_new(05) = 305;       rd_new(05) = 0;     rt_new(05) = 1-rt(05);  % NOTE SIGN CHANGE!
rx_new(06) = 192;               ry_new(06) = 305;       rd_new(06) = 0;     rt_new(06) = rt(06);
rx_new(07) = 188;               ry_new(07) = 279;       rd_new(07) = 0;     rt_new(07) = 1-rt(07);  % NOTE SIGN CHANGE!
rx_new(08) = 188;               ry_new(08) = 278;       rd_new(08) = 0;     rt_new(08) = 1-rt(08);  % NOTE SIGN CHANGE!

rx_new(09) = 172;               ry_new(09) = 236;
rx_new(10) = 172;               ry_new(10) = 237;
rx_new(21) = 192;               ry_new(21) = 298;      
rx_new(22) = 192;               ry_new(22) = 304;
rx_new(23) = 194;               ry_new(23) = 274;
rx_new(24) = 46;                ry_new(24) = 136;
rx_new(25) = 45;                ry_new(25) = 137;
rx_new(26) = 45;                ry_new(26) = 138;
rx_new(27) = 45;                ry_new(27) = 139;
rx_new(28:32) = 10;             ry_new(28:32) = 111:115;

% % Show u rivers
% i = find(rd_new==0);
% figure; 
%   subplot(1,2,1); hold on;
%     imagesc(xr-0.5,yr-0.5,mr'); set(gca,'ydir','normal'); axis tight; box on;
%     plot(rx(i)-0.5,ry(i),'or');
%     quiver(rx(i)-0.5,ry(i),sign(rt(i)),zeros(numel(i),1),1e-2,'-k');
%     for j=i; text(rx(j)-0.5,ry(j)+0.25,num2str(j)); end
%   hold off;
%   subplot(1,2,2); hold on;
%     imagesc(xr_new-0.5,yr_new-0.5,mr_new'); set(gca,'ydir','normal'); axis tight; box on;
%     plot(rx_new(i)-0.5,ry_new(i),'or');
%     quiver(rx_new(i)-0.5,ry_new(i),sign(rt_new(i)),zeros(numel(i),1),1e-2,'-k');
%     for j=i; text(rx_new(j)-0.5,ry_new(j)+0.25,num2str(j)); end
%   hold off;
  
% Fixes for y
rx_new(01) = 190;       ry_new(01) = 307;
rx_new(02) = 191;       ry_new(02) = 307;
rx_new(11) = 165;       ry_new(11) = 243;
rx_new(12) = 166;       ry_new(12) = 243;
rx_new(13) = 165;       ry_new(13) = 251;
rx_new(14) = 166;       ry_new(14) = 251;
rx_new(15) = 203;       ry_new(15) = 250;
rx_new(16) = 204;       ry_new(16) = 250;
rx_new(17) = 232;       ry_new(17) = 242;
rx_new(18) = 233;       ry_new(18) = 242;
rx_new(19) = 080;       ry_new(19) = 140;
rx_new(20) = 081;       ry_new(20) = 140;

% Show v rivers
i = find(rd_new==1);
figure('units','normalized','outerposition',[0 0 1 1]);
  subplot(1,2,1); hold on;
    imagesc(xr-0.5,yr-0.5,mr'); set(gca,'ydir','normal'); axis tight; box on;
    plot(rx(i),ry(i)-0.5,'or');
    quiver(rx(i),ry(i)-0.5,zeros(numel(i),1),sign(rt(i)),1e-2,'-k');
    for j=i; text(rx(j),ry(j)-0.25,num2str(j)); end
  hold off;
  subplot(1,2,2); hold on;
    imagesc(xr_new-0.5,yr_new-0.5,mr_new'); set(gca,'ydir','normal'); axis tight; box on;
    plot(rx_new(i),ry_new(i)-0.5,'or');
    quiver(rx_new(i),ry_new(i)-0.5,zeros(numel(i),1),sign(rt_new(i)),1e-2,'-k');
    for j=i; text(rx_new(j),ry_new(j)-0.25,num2str(j)); end
  hold off;
  
% % Generate new river file
% newFile = 'test_rivers_newgrid.nc';
% info=ncinfo(rvrFile);
% info.Dimensions(1).Length = 330;
% info.Dimensions(2).Length = 329;
% info.Dimensions(3).Length = 330;
% info.Dimensions(4).Length = 310;
% info.Dimensions(5).Length = 310;
% info.Dimensions(6).Length = 309;
% ncwriteschema(newFile,info);
% for i=1:numel(info.Variables);  ncwrite(newFile,info.Variables(i).Name,ncread(rvrFile,info.Variables(i).Name)); end

% Once satisfied with new positions, write them to new file
ncwrite(newFile,'river_Xposition',rx_new);
ncwrite(newFile,'river_Eposition',ry_new);
ncwrite(newFile,'river_direction',rd_new);
ncwrite(newFile,'river_transport',repmat(rt_new(:),[1 info.Dimensions(end).Length]).*ncread(rvrFile,'river_transport'));

