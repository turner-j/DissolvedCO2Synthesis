% Clear workspace
close all;clear; clc;

USALQ = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','US-ALQ_DT');
USLos = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','US-Los_DT');

run ImportUSEvMREddyProcOut.m
USGCE = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','US-GCE_DT');

MBPPW1 = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','MBPPW1FCO2_interp');
MBPPW2 = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','MBPPW2FCO2_interp');

addpath 'E:\SanDiskSecureAccess\all MATLAB files\Residence paper\USEvM'

USICs = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','USICs_DT');
USICs2 = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','USICs_DTinterp');

USOWC = readtable("E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls",'FileType','spreadsheet','Sheet','US-OWC_DT');
USOWC2 = readtable("E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls",'FileType','spreadsheet','Sheet','US-OWC_DT_interp');
USOWC = USOWC(6:end,:);
OWCsoilCO2 = USOWC(:,12:end);
USOWC.pCO2 = mean(table2array(OWCsoilCO2),2,'omitnan');% take the mean of each row

APEX = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','APEX');% US-BZF

USEDN = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','US-EDN_DT');

load USMyb.mat

%% Eliminate unrealistic values
USLos.pCO2(USLos.pCO2<413)=NaN;
USALQ.pCO2(USALQ.pCO2<415)=NaN;
USMyb2.pCO2(USMyb2.pCO2<405.0)=NaN;
USMybREddyProcout.pCO2(USMybREddyProcout.pCO2<405.0)=NaN;
USMyb2 = USMyb2(145:end,:);

%% Prairie Potholes and Karsts
t = tiledlayout(2,2);
title(t,'Seasonal R_{eco} Cycling','FontSize',20)

nexttile(3)
title('Prairie Potholes and Karsts')
hold on
monthlyReco = accumarray(month(MBPPW1.DATE),MBPPW1.Reco,[],@nanmean,NaN);
[~,id] = unique(month(MBPPW1.DATE));
plot(day(MBPPW1.DATE(id),'dayofyear'),monthlyReco(~isnan(monthlyReco)),'o-','Color',[0.0588 0.4902 0.5216],'MarkerFaceColor',[0.0588 0.4902 0.5216],'MarkerEdgeColor',[0.0588 0.4902 0.5216],'LineWidth',1.5,'MarkerSize',6)

monthlyReco = accumarray(month(MBPPW2.DATE),MBPPW2.Reco,[],@nanmean,NaN);
[~,id] = unique(month(MBPPW2.DATE));
plot(day(MBPPW2.DATE(id),'dayofyear'),monthlyReco(~isnan(monthlyReco)),'o-','Color',[0.5961 0.8706 0.8902],'MarkerFaceColor',[0.5961 0.8706 0.8902],'MarkerEdgeColor',[0.5961 0.8706 0.8902],'LineWidth',1.5,'MarkerSize',6)

plot(1,nan,'o-','Color',[0.5608 0.6588 0.5843],'MarkerFaceColor',[0.5608 0.6588 0.5843],'MarkerEdgeColor','k','LineWidth',1.5,'MarkerSize',6)% US-ALQ
plot(1,nan,'o-','Color',[0.0588 0.3216 0.1294],'MarkerFaceColor',[0.0588 0.3216 0.1294],'MarkerEdgeColor',[0.0588 0.3216 0.1294],'LineWidth',1.5,'MarkerSize',6)% SE-Deg
plot(1,nan,'o--k','MarkerFaceColor',[0.2039 0.6196 0.3137],'MarkerEdgeColor','k','LineWidth',1.5,'MarkerSize',5)% CA-DBB
plot(1,nan,'o-','Color',[0.3686 1.0000 0.4941],'MarkerFaceColor',[0.3686 1.0000 0.4941],'MarkerEdgeColor',[0.3686 1.0000 0.4941],'LineWidth',1.5,'MarkerSize',6)% Troll
plot(1,nan,'o-','Color',[0.6706 0.7882 0.7176],'MarkerFaceColor',[0.6706 0.7882 0.7176],'MarkerEdgeColor',[0.6706 0.7882 0.7176],'LineWidth',1.5,'MarkerSize',6)% US-Myb
plot(1,nan,'o-','Color',[0.6627 0.9882 0.0118],'MarkerFaceColor',[0.6627 0.9882 0.0118],'MarkerEdgeColor','k','LineWidth',1.5,'MarkerSize',5) % US-Los (PW)
plot(1,nan,'o-','Color',[0.3333 0.4902 0.0235],'MarkerFaceColor',[0.3333 0.4902 0.0235],'MarkerEdgeColor',[0.3333 0.4902 0.0235],'LineWidth',1.5,'MarkerSize',6)% US-ALQ (PW)
plot(1,nan,'o-k','MarkerFaceColor',[0.7373 0.8196 0.2118],'MarkerEdgeColor',[0.7373 0.8196 0.2118],'LineWidth',1.5,'MarkerSize',5)% SE-Deg (PW)

ylabel('R_{eco} (\mumol CO_2 m^-^2s^-^1)')
hold off
xlim([0 370])
xline(182);
xlabel('Day of Year')
set(gca,'FontSize',20)
box on
leg = legend({'MBPPW1','MBPPW2','US-ALQ','SE-Deg','CA-DBB','Troll','US-Myb','US-Los (PW)','US-ALQ (PW)','SE-Deg (PW)'},...
    'FontSize',20,'NumColumns',4);
leg.Layout.Tile = 1;

%% fens, bogs, and marshes
load ALQStream2021.mat
CADBB = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','CA-DBB');% not the interpolated data
SEDeg = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','SEDeg_DT');
PorewaterSEDeg = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','PorewaterSEDeg_DT');
Troll = readtable("E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls",'FileType','spreadsheet','Sheet','Troll');
USMyb2 = USMyb2(145:end,:);

ax2 = nexttile(2);
hold on
monthlyReco = accumarray(month(X.TIMESTAMP_END),X.RECO_F,[],@nanmean,NaN);
[~,id] = unique(month(X.TIMESTAMP_END));
plot(day(X.TIMESTAMP_END(id),'dayofyear'),monthlyReco(~isnan(monthlyReco)),'o-','Color',[0.5608 0.6588 0.5843],'MarkerFaceColor',[0.5608 0.6588 0.5843],'MarkerEdgeColor','k','LineWidth',1.5,'MarkerSize',6)

monthlyReco = accumarray(month(SEDeg.DateTime),SEDeg.Reco_DT_uStar,[],@nanmean,NaN);
[~,id] = unique(month(SEDeg.DateTime));
plot(day(SEDeg.DateTime(id),'dayofyear'),monthlyReco(~isnan(monthlyReco)),'o-','Color',[0.0588 0.3216 0.1294],'MarkerFaceColor',[0.0588 0.3216 0.1294],'MarkerEdgeColor',[0.0588 0.3216 0.1294],'LineWidth',1.5,'MarkerSize',6)

monthlyReco = accumarray(month(CADBB.DATE),CADBB.Reco,[],@nanmean,NaN);
[~,id] = unique(month(CADBB.DATE));
plot(day(CADBB.DATE(id),'dayofyear'),monthlyReco(~isnan(monthlyReco)),'o--k','MarkerFaceColor',[0.2039 0.6196 0.3137],'MarkerEdgeColor','k','LineWidth',1.5,'MarkerSize',6)

monthlyReco = accumarray(month(Troll.DateTime),Troll.Reco_uStar,[],@nanmean,NaN);
[~,id] = unique(month(Troll.DateTime));
plot(day(Troll.DateTime(id),'dayofyear'),monthlyReco(~isnan(monthlyReco)),'o-','Color',[0.3686 1.0000 0.4941],'MarkerFaceColor',[0.3686 1.0000 0.4941],'MarkerEdgeColor',[0.3686 1.0000 0.4941],'LineWidth',1.5,'MarkerSize',6)

monthlyReco = accumarray(month(USMyb2.DateTime),USMyb2.Reco_DT_uStar,[],@nanmean,NaN);
[~,id] = unique(month(USMyb2.DateTime));
plot(day(USMyb2.DateTime(id),'dayofyear'),monthlyReco,'o-','Color',[0.6706 0.7882 0.7176],'MarkerFaceColor',[0.6706 0.7882 0.7176],'MarkerEdgeColor',[0.6706 0.7882 0.7176],'LineWidth',1.5,'MarkerSize',6)

xline(201.4);
xticklabels({''})
xticks(0:50:350)
hold off
xlim([0 370])
box on
title('Fens, Bogs, and Marshes')
set(gca,'FontSize',20)
ax2.YAxisLocation = 'right';
%% Porewater
ax = nexttile(4);
title('Porewater')
hold on
monthlyReco = accumarray(month(USLos.TIMESTAMP_END),USLos.Reco,[],@nanmean,NaN);
[~,id] = unique(month(USLos.TIMESTAMP_END));
plot(day(USLos.TIMESTAMP_END(id),'dayofyear'),monthlyReco(~isnan(monthlyReco)),'o-','Color',[0.6627 0.9882 0.0118],'MarkerFaceColor',[0.6627 0.9882 0.0118],'MarkerEdgeColor','k','LineWidth',1.5,'MarkerSize',6)

monthlyReco = accumarray(month(USALQ.TIMESTAMP_END),USALQ.RECO_F,[],@nanmean,NaN);
[~,id] = unique(month(USALQ.TIMESTAMP_END));
plot(day(USALQ.TIMESTAMP_END(id),'dayofyear'),monthlyReco(~isnan(monthlyReco)),'o-','Color',[0.3333 0.4902 0.0235],'MarkerFaceColor',[0.3333 0.4902 0.0235],'MarkerEdgeColor',[0.3333 0.4902 0.0235],'LineWidth',1.5,'MarkerSize',6)

monthlyReco = accumarray(month(PorewaterSEDeg.DateTime),PorewaterSEDeg.Reco_DT_uStar,[],@nanmean,NaN);
[~,id] = unique(month(PorewaterSEDeg.DateTime));
plot(day(PorewaterSEDeg.DateTime(id),'dayofyear'),monthlyReco(~isnan(monthlyReco)),'o-k','MarkerFaceColor',[0.7373 0.8196 0.2118],'MarkerEdgeColor',[0.7373 0.8196 0.2118],'LineWidth',1.5,'MarkerSize',6)

box on
ax.YAxisLocation = 'right';
t.TileSpacing = 'compact';
t.Padding = 'compact';
xline(182);
xlabel('Day of Year')
xlim([0 370])
set(gca,'FontSize',20)

