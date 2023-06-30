clear; clc;
run('E:\SanDiskSecureAccess\all MATLAB files\Residence paper\colorpalette.m\colorpalette.m')
addpath('E:\SanDiskSecureAccess\all MATLAB files\allequash')
USGCE = readtable("E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls",'FileType','spreadsheet','Sheet','US-GCEinterp');
tides = readtable('E:\SanDiskSecureAccess\all MATLAB files\Sapelo\Sapelo SST\saphdwq2021.csv');% Hunt Dock

%% Plotting
xx = tides.DateTimeStamp;
yy = tides.cDepth;
tide = interp1(xx,yy,USGCE.TIMESTAMP_END);% Interpolate tide data to match co2lamp times

h = figure();

right_color = [0 0 0];
left_color = [0 0 0];
set(h,'defaultAxesColorOrder',[left_color; right_color]);
c = tide(day(USGCE.TIMESTAMP_END,'dayofyear')==222);

yyaxis right
hold on 
plot(USGCE.TIMESTAMP_END(day(USGCE.TIMESTAMP_END,'dayofyear')==222),tide(day(USGCE.TIMESTAMP_END,'dayofyear')==222),'-k','LineWidth',2)
hold off
ylabel('Tide Ht (m)')
yyaxis left
plot(USGCE.TIMESTAMP_END(day(USGCE.TIMESTAMP_END,'dayofyear')==222),USGCE.pCO2(day(USGCE.TIMESTAMP_END,'dayofyear')==222),'dk','MarkerFaceColor',[0.8471 0.8235 0.9686],'MarkerEdgeColor','k','MarkerSize',8)
ylabel('CO_2 {\itaq} (ppmv)')

xlabel('Hour')
legend('P_{CO2}','Tide Ht')
set(gca,'FontSize',12)
