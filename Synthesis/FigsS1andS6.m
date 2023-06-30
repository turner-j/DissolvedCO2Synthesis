% Clear workspace
close all;clear; clc;
addpath(genpath(pwd))

load("USALQSWPCO2.mat")
load USEDN.mat
load('LochValepCO2.mat')
LochValeCO2 = Mastetal1998LochValeCO2flux.SaturatedAverageDailyCO2Fluxmmolm2d1;% Saturated soil FCO2
LochValeCO2 = LochValeCO2.*1000./86400;% Convert from units per day to seconds
load("LVDailyData.mat");LVDailyData = LVDailyData(LVDailyData.Site_ID==401723105400000,:);
run ImportLVMetData.m;LochValeClimate = WY1992to2019LochValeClimate((WY1992to2019LochValeClimate.station_name=='Andrews Creek weather station')&(WY1992to2019LochValeClimate.measurement_height==4),:);

USLos2 = readtable("E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls",'FileType','spreadsheet','Sheet','US-Los');
USALQ2 = readtable("E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls",'FileType','spreadsheet','Sheet','US-ALQ_DT');

run ImportBigCypresspCO2.m

PLMdiel = readtable('E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls','Sheet','PLM');
PHMdiel = readtable('E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls','Sheet','PHM');

USHB1diel = readtable('E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls','Sheet','USHB1');

PorewaterSEDeg = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','PorewaterSEDeg_DT');
load USMyb.mat;USMyb2 = USMyb2(145:end,:);

USGCE2 = readtable("E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls",'FileType','spreadsheet','Sheet','US-GCEinterp');
USEvM2 = readtable("E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls",'FileType','spreadsheet','Sheet','USEvM_DT');

%% Calculate monthly average snow depth
ind = find(ismember(LochValeClimate.timestamp,LVDailyData.Date));
LochValeClimate = LochValeClimate(ind(1):ind(end),:);

s = string(LochValeClimate.Snow_depth);LochValeClimate.Snow_depth = double(s);
monthlysnow = accumarray(month(LochValeClimate.timestamp),LochValeClimate.Snow_depth,[],@(x)mean(x,'omitnan'),NaN);
%% Eliminate unrealistic values
USLos2.pCO2(USLos2.pCO2<413)=NaN;
USALQ2.pCO2(USALQ2.pCO2<415)=NaN;
USMyb2.pCO2(USMyb2.pCO2<405.0)=NaN;
USMyb2 = USMyb2(145:end,:);
%% Plot Tidal Site pCO2 & NEE
t = tiledlayout(2,2);
title(t,'Diel P_{CO2} Cycling','FontSize',20)

% Alpine Sites
LVDiel = readtable("E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls",'FileType','spreadsheet','Sheet','LV_DT');
LVDiel.Snow_depth = str2double(LVDiel.Snow_depth);
ix = find(ismember(LVDiel.timestamp,LochValeClimate.timestamp));
LVDiel.Snow_depth = LochValeClimate.Snow_depth(ix,:);

nexttile(1)
dailypCO2 = accumarray(hour(LVDiel.timestamp)+1,LVDiel.pCO2,[],@(x)mean(x,'omitnan'));
plot(0:23,dailypCO2,'o','MarkerFaceColor',[0.2118 0.3529 0.4784],'MarkerEdgeColor',[0.2118 0.3529 0.4784],'MarkerSize',8)
box on
xlim([0 24])
title('Alpine')
set(gca,'FontSize',20)
ylabel('CO_2({\itaq}) (ppmv)')
xticklabels({''})
xticks(0:5:20)

Daily = dailypCO2;

% fens, bogs, and marshes
run ImportTrollFCO2.m
run ImportTrollpCO2.m

ax3 = nexttile(2);
hold on

YKDB2 = readtable("E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls",'FileType','spreadsheet','Sheet','YKDB');
YKDUB2 = readtable("E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls",'FileType','spreadsheet','Sheet','YKDUB');
[hourlypCO2,~] = dielNEE(YKDUB2.datetime,YKDUB2.pCO2);

[hrlypCO2,~,~] = groupsummary(YKDUB2.pCO2,round(hour(YKDUB2.datetime))+1,'mean');
plot(0:.5:23.5,hourlypCO2,'o','MarkerFaceColor',[0.0353 0.9098 0.2706],'MarkerEdgeColor','k','MarkerSize',8)

[hourlypCO2,~] = dielNEE(USALQ.TIMESTAMP_END,USALQ.SWpCO2);
[dailypCO2,~,~] = groupsummary(USALQ.SWpCO2,round(hour(USALQ.TIMESTAMP_END))+1,'mean');

plot(0:.5:23.5,hourlypCO2,'s','MarkerFaceColor',[0.5608 0.6588 0.5843],'MarkerEdgeColor','k','MarkerSize',8)

Daily = [Daily,hrlypCO2,dailypCO2]; 

[hourlypCO2,~] = dielNEE(USMybREddyProcout.DateTime,USMybREddyProcout.pCO2);
plot(0:.5:23.5,hourlypCO2,'o-k','MarkerFaceColor',[0.6706 0.7882 0.7176],'MarkerEdgeColor','k','MarkerSize',8)

[hrlypCO2,~,~] = groupsummary(USMybREddyProcout.pCO2,round(hour(USMybREddyProcout.DateTime))+1,'mean');

[hourlypCO2,~] = dielNEE(YKDB2.datetime,YKDB2.pCO2);
plot(0:.5:23.5,hourlypCO2,'o','MarkerFaceColor',[0.0157 0.6902 0.1961],'MarkerEdgeColor',[0.0157 0.6902 0.1961],'MarkerSize',8)

[dailypCO2,~,~] = groupsummary(YKDB2.pCO2,round(hour(YKDB2.datetime))+1,'mean');

Daily = [Daily,hrlypCO2,dailypCO2];

run ImportImnavaitpCO2.m
[meanpCO2,~,BC] = groupsummary(ImnavaitCreekpCO2.PCO2_uatm(~isnat(ImnavaitCreekpCO2.Time_hr_dst)),hour(ImnavaitCreekpCO2.Time_hr_dst(~isnat(ImnavaitCreekpCO2.Time_hr_dst))),'mean');
timestamps = unique(hour(ImnavaitCreekpCO2.Time_hr_dst(~isnat(ImnavaitCreekpCO2.Time_hr_dst))));

if (sum(BC>=5)>=10)
    plot(timestamps(~isnan(meanpCO2)),meanpCO2(~isnan(meanpCO2)),'o','MarkerFaceColor',[0.7569 0.9686 0.8118],'MarkerEdgeColor','k','MarkerSize',8)
end

Daily = [Daily,meanpCO2];

hold off
xticklabels({''})
xticks(0:5:20)
box on
ax3.YAxisLocation = 'right';
xlim([0 24])
title('Fens, Bogs, and Marshes')
set(gca,'FontSize',20)
addpath('E:\SanDiskSecureAccess\all MATLAB files\Residence paper\ImnavaitCreek\')
run ImportImnavaitREddyOut.m

% Porewater
load SEDegPorewater.mat
SEDegREddyProcOut2 = SEDegREddyProcOut2(7915:77256,:);

USLos2 = USLos2((ismember(day(USLos2.TIMESTAMP_END,'dayofyear'),229:366)&ismember(year(USLos2.TIMESTAMP_END),2020)),:);
USLos2.pCO2(USLos2.pCO2<413)=NaN;
USALQ2.pCO2(USALQ2.pCO2<415)=NaN;

ax2 = nexttile(4);
yyaxis left
hold on
[hourlypCO2,~] = dielNEE(USLos2.TIMESTAMP_END,USLos2.pCO2);
plot(0:.5:23.5,hourlypCO2,'o','MarkerFaceColor',[0.6627 0.9882 0.0118],'MarkerEdgeColor','k','MarkerSize',8)%US-Los (PW)
[hourlypCO2,~] = dielNEE(USALQ2.TIMESTAMP_END,USALQ2.pCO2);
plot(0:.5:23.5,hourlypCO2,'o','MarkerFaceColor',[0.3333 0.4902 0.0235],'MarkerEdgeColor',[0.3333 0.4902 0.0235],'MarkerSize',8)% US-ALQ (PW)
hold off
ax2.XColor = [0 0 0];
ax2.YColor = [0 0 0];
xlim([0 24])
xlabel('Hour')

yyaxis right
hold on
[hourlypCO2,~] = dielNEE(SEDegREddyProcOut2.DateTime,SEDegREddyProcOut2.pCO2);
plot(0:.5:23.5,hourlypCO2,'o','MarkerFaceColor',[0.7373 0.8196 0.2118],'MarkerEdgeColor',[0.7373 0.8196 0.2118],'MarkerSize',8)%SE-Deg (PW)

% Use blanks for the sites without data in order to autofill the legend
plot(NaN,NaN,'o','MarkerFaceColor',[0.2118 0.3529 0.4784],'MarkerEdgeColor',[0.2118 0.3529 0.4784],'MarkerSize',8)%LV
plot(NaN,NaN,'o','MarkerFaceColor',[0.0157 0.6902 0.1961],'MarkerEdgeColor',[0.0157 0.6902 0.1961],'MarkerSize',8)%YKDB
plot(NaN,NaN,'o','MarkerFaceColor',[0.0353 0.9098 0.2706],'MarkerEdgeColor','k','MarkerSize',8)% YKDUB
plot(NaN,NaN,'s','MarkerFaceColor',[0.5608 0.6588 0.5843],'MarkerEdgeColor','k','MarkerSize',8)% US-ALQ Stream
plot(NaN,NaN,'o','MarkerFaceColor',[0.7569 0.9686 0.8118],'MarkerEdgeColor','k','MarkerSize',8)%US-ICs
plot(NaN,NaN,'o-k','MarkerFaceColor',[0.6706 0.7882 0.7176],'MarkerEdgeColor','k','MarkerSize',8)%USMyb
hold off
box on

ax2.XColor = [0 0 0];
ax2.YColor = [0 0 0];
xlim([0 24])
title('Porewater')
set(gca,'FontSize',20)

leg = legend({'US-Los (PW)','US-ALQ (PW)','SE-Deg (PW)','LV','YKD 1','YKD 2','US-ALQ','US-ICs','US-Myb'},...
    'FontSize',20,'NumColumns',4);
leg.Layout.Tile = 3;

t.TileSpacing = 'compact';
t.Padding = 'compact';
%% Tidal wetland subplots
figure()

z = tiledlayout(4,1);
title(z,'Tidal P_{CO2} Timeseries','FontSize',20)
ylabel(z,'CO_2({\itaq}) (ppmv)','FontSize',20)

nexttile(1)
hold on 
USGCE2 = USGCE2((ismember(day(USGCE2.TIMESTAMP_END,'dayofyear'),203:209)),:);
plot(USGCE2.TIMESTAMP_END(~isnan(USGCE2.pCO2)),USGCE2.pCO2(~isnan(USGCE2.pCO2)),'o-k','MarkerFaceColor',[0.8471 0.8235 0.9686],'MarkerEdgeColor','k','MarkerSize',8)

plot(USGCE2.TIMESTAMP_END(1),NaN,'o','MarkerFaceColor',[0.1373 0.0863 0.4784],'MarkerEdgeColor','k','MarkerSize',8)%US-EvM
plot(USGCE2.TIMESTAMP_END(1),NaN,'o-','Color',[0.7804 0.4314 1.0000],'LineWidth',1,'MarkerFaceColor',[0.7804 0.4314 1.0000],'MarkerEdgeColor','k','MarkerSize',8)%US-HB1
plot(USEDN.DateTime(1),NaN,'^-','Color',[0.7882 0.6000 1.0000],'LineWidth',1,'MarkerFaceColor',[0.7882 0.6000 1.0000],'MarkerEdgeColor','k','MarkerSize',8)

hold off
box on
set(gca,'FontSize',20)
leg = legend({'US-GCE','US-EvM','US-HB1','US-EDN'},'Location','northwest','FontSize',15);

nexttile(2)
plot(USEvM2.DateTime,USEvM2.pCO2,'o','MarkerFaceColor',[0.1373 0.0863 0.4784],'MarkerEdgeColor','k','MarkerSize',8)
set(gca,'FontSize',20)

nexttile(3)
ix = find(ismember(day(USHB1diel.Time,'dayofyear'),225:231));
plot(USHB1diel.Time(ix),USHB1diel.pCO2(ix),'o-','Color',[0.7804 0.4314 1.0000],'LineWidth',1,'MarkerFaceColor',[0.7804 0.4314 1.0000],'MarkerEdgeColor','k','MarkerSize',8)
set(gca,'FontSize',20)

nexttile(4)
ix = find(ismember(day(USEDN.DateTime,'dayofyear'),121:127));
plot(USEDN.DateTime(ix),USEDN.CO2(ix),'^-','Color',[0.7882 0.6000 1.0000],'LineWidth',1,'MarkerFaceColor',[0.7882 0.6000 1.0000],'MarkerEdgeColor','k','MarkerSize',8)
set(gca,'FontSize',20)
z.TileSpacing = 'compact';
z.Padding = 'compact';
