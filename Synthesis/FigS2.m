% Clear workspace
close all;clear; clc;
run ImportUSEvMREddyProcOut.m
load("USALQSWPCO2.mat")

addpath 'E:\SanDiskSecureAccess\all MATLAB files\Residence paper\USEvM'
addpath('E:\SanDiskSecureAccess\all MATLAB files\Residence paper\YKD')

run ImportBigCypresspCO2.m

PLMdiel = readtable('E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls','Sheet','PLM');
PHMdiel = readtable('E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls','Sheet','PHM');

USHB1diel = readtable('E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls','Sheet','USHB1');

load USMyb.mat;USMyb2 = USMyb2(145:end,:);
load USEDN.mat

USLos2 = readtable("E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls",'FileType','spreadsheet','Sheet','US-Los');
USALQ2 = readtable("E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls",'FileType','spreadsheet','Sheet','US-ALQ_DT');

load SEDegREddyProcOut2.mat

run ImportTrollFCO2.m
run ImportTrollpCO2.m

addpath('E:\SanDiskSecureAccess\all MATLAB files\allequash')
USGCE2 = readtable("E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls",'FileType','spreadsheet','Sheet','US-GCEinterp');
USEvM2 = readtable("E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls",'FileType','spreadsheet','Sheet','USEvM_DT');

addpath('E:\SanDiskSecureAccess\all MATLAB files\Residence paper\ImnavaitCreek\')
run ImportImnavaitREddyOut.m
%% Eliminate unrealistic values
USALQ.SWpCO2(USALQ.SWpCO2<415)=NaN;
USMybREddyProcout.pCO2(USMybREddyProcout.pCO2<405.0)=NaN;
%% Plot Reco
t = tiledlayout(2,2);
title(t,'Diel R_{eco} Cycling','FontSize',20)

% fens, bogs, marshes
ax3 = nexttile(2);
hold on

for j= 2:width(USICsREddyProcOut)
    USICsREddyProcOut.(j)(USICsREddyProcOut.(j)<=-9999) = NaN;
end

[meanReco,~,BC] = groupsummary(USICsREddyProcOut.Reco_DT_uStar(~isnat(USICsREddyProcOut.DateTime)),hour(USICsREddyProcOut.DateTime(~isnat(USICsREddyProcOut.DateTime))),'mean');

if (sum(BC>=5)>=10)
    plot(unique(hour(USICsREddyProcOut.DateTime(~isnat(USICsREddyProcOut.DateTime)))),meanReco,'ok','MarkerFaceColor',[0.7569 0.9686 0.8118],'MarkerEdgeColor','k','MarkerSize',8)
end

[hourlyReco,~] = dielNEE(TrollFCO2.DateTime,TrollFCO2.Reco_uStar);
plot(0:.5:23.5,hourlyReco,'o','MarkerFaceColor',[0.3686 1.0000 0.4941],'MarkerEdgeColor',[0.3686 1.0000 0.4941],'MarkerSize',8)

USMybREddyProcout = USMybREddyProcout(6938:end,:);
[hourlyReco,~] = dielNEE(USMybREddyProcout.DateTime,USMybREddyProcout.Reco_DT_uStar);
plot(0:.5:23.5,hourlyReco,'o','MarkerFaceColor',[0.6706 0.7882 0.7176],'MarkerEdgeColor',[0.6706 0.7882 0.7176],'MarkerSize',8)

hold off
ax3.YAxisLocation = 'right';
xticklabels({''})
xticks(0:5:20)
box on
xlim([0 24])
title('Fens, Bogs, and Marshes')
set(gca,'FontSize',20)
addpath('E:\SanDiskSecureAccess\all MATLAB files\Residence paper\ImnavaitCreek\')
run ImportImnavaitREddyOut.m

% Tidal
nexttile(1);
title('Tidal')

hold on
[hourlyReco,~] = dielNEE(USGCE2.TIMESTAMP_END,USGCE2.Reco);
plot(0:.5:23.5,hourlyReco,'s','MarkerFaceColor',[0.8471 0.8235 0.9686],'MarkerEdgeColor','k','MarkerSize',8)
[hourlyReco,~] = dielNEE(USEvMREddyProcOut.DateTime,USEvMREddyProcOut.Reco_DT_uStar);
plot(0:.5:23.5,hourlyReco,'o','MarkerFaceColor',[0.1373 0.0863 0.4784],'MarkerEdgeColor','k','MarkerSize',8)
[hourlyReco,~] = dielNEE(PLMdiel.TIMESTAMP_END,PLMdiel.Reco_DT_uStar);
plot(0:.5:23.5,hourlyReco,'o','MarkerFaceColor',[0.5020 0.3608 0.5882],'MarkerEdgeColor',[0.5020 0.3608 0.5882],'MarkerSize',8)
[hourlyReco,~] = dielNEE(USHB1diel.Time,USHB1diel.Reco);
plot(0:.5:23.5,hourlyReco,'o','MarkerFaceColor',[0.7804 0.4314 1.0000],'MarkerEdgeColor',[0.7804 0.4314 1.0000],'MarkerSize',8)
[hourlyReco,~] = dielNEE(USEDN.DateTime,USEDN.Reco_DT_uStar);
plot(0:.5:23.5,hourlyReco,'^','MarkerFaceColor',[0.7882 0.6000 1.0000],'MarkerEdgeColor','k','MarkerSize',8)

% Use blanks for the sites without data in order to autofill the legend
plot(NaN,NaN,'o','MarkerFaceColor',[0.7569 0.9686 0.8118],'MarkerEdgeColor','k','MarkerSize',8)%US-ICs
plot(NaN,NaN,'o','MarkerFaceColor',[0.6706 0.7882 0.7176],'MarkerEdgeColor',[0.6706 0.7882 0.7176],'MarkerSize',8)%USMyb
plot(NaN,NaN,'o','MarkerFaceColor',[0.3686 1.0000 0.4941],'MarkerEdgeColor',[0.3686 1.0000 0.4941],'MarkerSize',8)% Troll.
plot(NaN,NaN,'o','MarkerFaceColor',[0.6627 0.9882 0.0118],'MarkerEdgeColor','k','MarkerSize',8)%US-Los (PW)
plot(NaN,NaN,'s','MarkerFaceColor',[0.3333 0.4902 0.0235],'MarkerEdgeColor',[0.3333 0.4902 0.0235],'MarkerSize',8)% US-ALQ (PW)
plot(NaN,NaN,'o','MarkerFaceColor',[0.7373 0.8196 0.2118],'MarkerEdgeColor',[0.7373 0.8196 0.2118],'MarkerSize',8)%SE-Deg (PW)

hold off
ylabel('R_{eco} (\mumol CO_2 m^-^2s^-^1)')
box on
xlim([0 24])
xlabel('Hour')
lg = legend({'US-GCE','US-EvM','PIE','US-HB1','US-EDN','US-ICs','US-Myb','Troll.','US-Los (PW)','US-ALQ (PW)','SE-Deg (PW)'},...
    'Location','southoutside','FontSize',20,'NumColumns',4);
set(gca,'FontSize',20)
lg.Layout.Tile = 3;

% Porewater
USLos2 = USLos2((ismember(day(USLos2.TIMESTAMP_END,'dayofyear'),229:366)&ismember(year(USLos2.TIMESTAMP_END),2020)),:);
USLos2.pCO2(USLos2.pCO2<413)=NaN;
USALQ2.pCO2(USALQ2.pCO2<415)=NaN;

ax = nexttile(4);
hold on
[hourlyReco,~] = dielNEE(USLos2.TIMESTAMP_END,USLos2.RECO_PI_F);
plot(0:.5:23.5,hourlyReco,'o','MarkerFaceColor',[0.6627 0.9882 0.0118],'MarkerEdgeColor','k','MarkerSize',8)
[hourlyReco,~] = dielNEE(USALQ2.TIMESTAMP_END,USALQ2.RECO_F);
plot(0:.5:23.5,hourlyReco,'s','MarkerFaceColor',[0.3333 0.4902 0.0235],'MarkerEdgeColor','none','MarkerSize',8)
[hourlyReco,~] = dielNEE(SEDegREddyProcOut2.DateTime,SEDegREddyProcOut2.Reco_DT_uStar);
plot(0:.5:23.5,hourlyReco,'o','MarkerFaceColor',[0.7373 0.8196 0.2118],'MarkerEdgeColor','none','MarkerSize',8)
hold off
xlim([0 24])
xlabel('Hour')
box on
ax.YAxisLocation = 'right';

title('Porewater')
set(gca,'FontSize',20)

t.TileSpacing = 'compact';
t.Padding = 'compact';

