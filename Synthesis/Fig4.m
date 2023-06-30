% Clear workspace
close all;clear; clc;
addpath 'E:\SanDiskSecureAccess\all MATLAB files\Residence paper\US-EDN'
addpath 'E:\SanDiskSecureAccess\all MATLAB files\Residence paper'\CADBB\
load USALQStreamDaily.mat
USALQ = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','US-ALQ_DT');
SEDeg = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','SEDeg_DT');
PorewaterSEDeg = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','PorewaterSEDeg_DT');

load USEDNPCO2_daily.mat
run ImportUSENDREddyProcOut.m

addpath 'E:\SanDiskSecureAccess\all MATLAB files\allequash'
for j= 2:width(USEDNREddyProcOut)
    USEDNREddyProcOut.(j)(USEDNREddyProcOut.(j)==-9999) = NaN;
    dailyave = dailyaverage2(USEDNREddyProcOut.(j));
    USEDN(:,j) = dailyave(:,1);
end

% add variable names, timestamp
USEDN = array2table(USEDN);
USEDN.Properties.VariableNames = USEDNREddyProcOut.Properties.VariableNames;
USEDN.DateTime = USEDNREddyProcOut.DateTime(dailyave(:,2).*48);

%% Eliminate unrealistic values
USALQ.pCO2(USALQ.pCO2<415)=NaN;
USALQ = USALQ(104:296,:);
%% Plotting
% PCO2
z = tiledlayout(3,3);

f = nexttile(1);
ylabel(z,'CO_2({\itaq}) (ppmv)','FontSize',20)

yyaxis left 
hold on
monthlyPCO2 = accumarray(month(SEDeg.DateTime),SEDeg.pCO2,[],@nanmean,NaN);
[~,id] = unique(month(SEDeg.DateTime));
[stderr,~,~] = groupsummary(SEDeg.pCO2,month(SEDeg.DateTime),'std');
plot(day(SEDeg.DateTime(id),'dayofyear'),monthlyPCO2(~isnan(monthlyPCO2)),'o-','Color',[0.7373 0.8196 0.2118],'MarkerFaceColor',[0.7373 0.8196 0.2118],'MarkerEdgeColor','k','LineWidth',2,'MarkerSize',8)
patch([day(SEDeg.DateTime(id),'dayofyear'); flipud(day(SEDeg.DateTime(id),'dayofyear'))], [monthlyPCO2(~isnan(monthlyPCO2))-stderr(:);  flipud(monthlyPCO2(~isnan(monthlyPCO2))+stderr(:))],[0.7373 0.8196 0.2118],'EdgeColor','none')
alpha(.3)
hold off
f.XColor = [0 0 0];
f.YColor = [0 0 0];

var(monthlyPCO2(~isnan(monthlyPCO2)))

yyaxis right
monthlyPCO2 = accumarray(month(PorewaterSEDeg.DateTime),PorewaterSEDeg.PorewaterCO2,[],@nanmean,NaN);
[~,id] = unique(month(PorewaterSEDeg.DateTime));
plot(day(PorewaterSEDeg.DateTime(id),'dayofyear'),monthlyPCO2(~isnan(monthlyPCO2)),'.--k','MarkerSize',20,'LineWidth',2)
legend('','','porewater','Location','northwest')

var(monthlyPCO2(~isnan(monthlyPCO2)))

f.XColor = [0 0 0];
f.YColor = [0 0 0];
xlim([0 370])
xticks(0:50:350)
xticklabels({''})
set(gca,'FontSize',20)
box on

g = nexttile(7);
yyaxis left
hold on
monthlyPCO2 = accumarray(month(USALQStreamDaily.Date_Time),USALQStreamDaily.pCO2_ppm,[],@nanmean,NaN);
[~,id] = unique(month(USALQStreamDaily.Date_Time));
[stderr,~,~] = groupsummary(USALQStreamDaily.pCO2_ppm,month(USALQStreamDaily.Date_Time),'std');
plot(day(USALQStreamDaily.Date_Time(id),'dayofyear'),monthlyPCO2(~isnan(monthlyPCO2)),'s-','Color',[0.3333 0.4902 0.0235],'MarkerFaceColor',[0.3333 0.4902 0.0235],'MarkerEdgeColor','k','LineWidth',2,'MarkerSize',10)
patch([day(USALQStreamDaily.Date_Time(id),'dayofyear'); flipud(day(USALQStreamDaily.Date_Time(id),'dayofyear'))], [monthlyPCO2(~isnan(monthlyPCO2))-stderr(:);  flipud(monthlyPCO2(~isnan(monthlyPCO2))+stderr(:))],[0.3333 0.4902 0.0235],'EdgeColor','none')
alpha(.3)

var(monthlyPCO2(~isnan(monthlyPCO2)))

hold off
g.XColor = [0 0 0];
g.YColor = [0 0 0];

yyaxis right
monthlyPCO2 = accumarray(month(USALQ.TIMESTAMP_END),USALQ.pCO2,[],@nanmean,NaN);
[~,id] = unique(month(USALQ.TIMESTAMP_END));
plot(day(USALQ.TIMESTAMP_END(id),'dayofyear'),monthlyPCO2(~isnan(monthlyPCO2)),'.--k','MarkerSize',20,'LineWidth',2)
xticks(0:50:350)
xlim([0 370])
set(gca,'FontSize',20)
box on
g.XColor = [0 0 0];
g.YColor = [0 0 0];
legend('','','porewater','Location','northwest')

var(monthlyPCO2(~isnan(monthlyPCO2)))

nexttile(4)
hold on
monthlyPCO2 = accumarray(month(USEDNPCO2DAY.timestamp_sel),USEDNPCO2DAY.CO2,[],@nanmean,NaN);
[~,id] = unique(month(USEDNPCO2DAY.timestamp_sel));
[stderr,~,~] = groupsummary(USEDNPCO2DAY.CO2,month(USEDNPCO2DAY.timestamp_sel),'std');
plot(day(USEDNPCO2DAY.timestamp_sel(id),'dayofyear'),monthlyPCO2(~isnan(monthlyPCO2)),'^-','Color',[0.7882 0.6000 1.0000],'MarkerFaceColor',[0.7882 0.6000 1.0000],'MarkerEdgeColor','k','LineWidth',2,'MarkerSize',10)
patch([day(USEDNPCO2DAY.timestamp_sel(id),'dayofyear'); flipud(day(USEDNPCO2DAY.timestamp_sel(id),'dayofyear'))], [monthlyPCO2(~isnan(monthlyPCO2))-stderr(:);  flipud(monthlyPCO2(~isnan(monthlyPCO2))+stderr(:))],[0.7882 0.6000 1.0000],'EdgeColor','none')
alpha(.3)
hold off
xticks(0:50:350)
xticklabels({''})
xlim([0 370])
set(gca,'FontSize',20)
box on

var(monthlyPCO2(~isnan(monthlyPCO2)))

%% GPP
ax = nexttile([3 1]);
hold on
monthlyGPP = accumarray(month(PorewaterSEDeg.DateTime),PorewaterSEDeg.GPP_DT_uStar,[],@nanmean,NaN);
[~,id] = unique(month(PorewaterSEDeg.DateTime));
plot(day(PorewaterSEDeg.DateTime(id),'dayofyear'),monthlyGPP(~isnan(monthlyGPP)),'o-','Color',[0.7373 0.8196 0.2118],'MarkerFaceColor',[0.7373 0.8196 0.2118],'MarkerEdgeColor','k','LineWidth',2,'MarkerSize',10)
[stderr,~,~] = groupsummary(PorewaterSEDeg.GPP_DT_uStar,month(PorewaterSEDeg.DateTime),'std');
patch([day(PorewaterSEDeg.DateTime(id),'dayofyear'); flipud(day(PorewaterSEDeg.DateTime(id),'dayofyear'))], [monthlyGPP(~isnan(monthlyGPP))-stderr(:);  flipud(monthlyGPP(~isnan(monthlyGPP))+stderr(:))],[0.7373 0.8196 0.2118],'EdgeColor','none')
alpha(.3)

var(monthlyGPP(~isnan(monthlyGPP)))

monthlyGPP = accumarray(month(USALQ.TIMESTAMP_END),USALQ.GPP_F,[],@nanmean,NaN);
[~,id] = unique(month(USALQ.TIMESTAMP_END));
plot(day(USALQ.TIMESTAMP_END(id),'dayofyear'),monthlyGPP(~isnan(monthlyGPP)),'s-','Color',[0.3333 0.4902 0.0235],'MarkerFaceColor',[0.3333 0.4902 0.0235],'MarkerEdgeColor','k','LineWidth',2,'MarkerSize',10)
[stderr,~,~] = groupsummary(USALQ.GPP_F,month(USALQ.TIMESTAMP_END),'std');
patch([day(USALQ.TIMESTAMP_END(id),'dayofyear'); flipud(day(USALQ.TIMESTAMP_END(id),'dayofyear'))], [monthlyGPP(~isnan(monthlyGPP))-stderr(:);  flipud(monthlyGPP(~isnan(monthlyGPP))+stderr(:))],[0.3333 0.4902 0.0235],'EdgeColor','none')
alpha(.3)

var(monthlyGPP(~isnan(monthlyGPP)))

monthlyGPP = accumarray(month(USEDN.DateTime),USEDN.GPP_DT_uStar,[],@nanmean,NaN);
[~,id] = unique(month(USEDN.DateTime));
plot(day(USEDN.DateTime(id),'dayofyear'),monthlyGPP(~isnan(monthlyGPP)),'^-','Color',[0.7882 0.6000 1.0000],'MarkerFaceColor',[0.7882 0.6000 1.0000],'MarkerEdgeColor','k','LineWidth',2,'MarkerSize',10)
[stderr,~,~] = groupsummary(USEDN.GPP_DT_uStar,month(USEDN.DateTime),'std');
patch([day(USEDN.DateTime(id),'dayofyear'); flipud(day(USEDN.DateTime(id),'dayofyear'))], [monthlyGPP(~isnan(monthlyGPP))-stderr(:);  flipud(monthlyGPP(~isnan(monthlyGPP))+stderr(:))],[0.7882 0.6000 1.0000],'EdgeColor','none')
alpha(.3)

var(monthlyGPP(~isnan(monthlyGPP)))

ylabel('GPP (\mumol CO_2 m^-^2s^-^1)')
hold off
xticks(0:50:350)
xlim([0 370])
set(gca,'FontSize',20)
box on
leg = legend({'SE-Deg (SW)','','US-ALQ (SW)','','US-EDN (SW)',''});
xlabel('Day of Year')
ax.YAxisLocation = 'right';

% Reco
ax2 = nexttile([3 1]);
hold on
monthlyReco = accumarray(month(PorewaterSEDeg.DateTime),PorewaterSEDeg.Reco_DT_uStar,[],@nanmean,NaN);
[~,id] = unique(month(PorewaterSEDeg.DateTime));
plot(day(PorewaterSEDeg.DateTime(id),'dayofyear'),monthlyReco(~isnan(monthlyReco)),'o-','Color',[0.7373 0.8196 0.2118],'MarkerFaceColor',[0.7373 0.8196 0.2118],'MarkerEdgeColor','k','LineWidth',2,'MarkerSize',10)
[stderr,~,~] = groupsummary(PorewaterSEDeg.Reco_DT_uStar,month(PorewaterSEDeg.DateTime),'std');
patch([day(PorewaterSEDeg.DateTime(id),'dayofyear'); flipud(day(PorewaterSEDeg.DateTime(id),'dayofyear'))], [monthlyReco(~isnan(monthlyReco))-stderr(:);  flipud(monthlyReco(~isnan(monthlyReco))+stderr(:))],[0.7373 0.8196 0.2118],'EdgeColor','none')
alpha(.3)

var(monthlyReco(~isnan(monthlyReco)))

monthlyReco = accumarray(month(USALQ.TIMESTAMP_END),USALQ.RECO_F,[],@nanmean,NaN);
[~,id] = unique(month(USALQ.TIMESTAMP_END));
plot(day(USALQ.TIMESTAMP_END(id),'dayofyear'),monthlyReco(~isnan(monthlyReco)),'s-','Color',[0.3333 0.4902 0.0235],'MarkerFaceColor',[0.3333 0.4902 0.0235],'MarkerEdgeColor','k','LineWidth',2,'MarkerSize',10)
[stderr,~,~] = groupsummary(USALQ.RECO_F,month(USALQ.TIMESTAMP_END),'std');
patch([day(USALQ.TIMESTAMP_END(id),'dayofyear'); flipud(day(USALQ.TIMESTAMP_END(id),'dayofyear'))], [monthlyReco(~isnan(monthlyReco))-stderr(:);  flipud(monthlyReco(~isnan(monthlyReco))+stderr(:))],[0.3333 0.4902 0.0235],'EdgeColor','none')
alpha(.3)

var(monthlyReco(~isnan(monthlyReco)))

monthlyReco = accumarray(month(USEDN.DateTime),USEDN.Reco_DT_uStar,[],@nanmean,NaN);
[~,id] = unique(month(USEDN.DateTime));
plot(day(USEDN.DateTime(id),'dayofyear'),monthlyReco(~isnan(monthlyReco)),'^-','Color',[0.7882 0.6000 1.0000],'MarkerFaceColor',[0.7882 0.6000 1.0000],'MarkerEdgeColor','k','LineWidth',2,'MarkerSize',10)
[stderr,~,~] = groupsummary(USEDN.Reco_DT_uStar,month(USEDN.DateTime),'std');
patch([day(USEDN.DateTime(id),'dayofyear'); flipud(day(USEDN.DateTime(id),'dayofyear'))], [monthlyReco(~isnan(monthlyReco))-stderr(:);  flipud(monthlyReco(~isnan(monthlyReco))+stderr(:))],[0.7882 0.6000 1.0000],'EdgeColor','none')
alpha(.3)

var(monthlyReco(~isnan(monthlyReco)))

xticks(0:50:350)
hold off
xlim([0 370])
box on
set(gca,'FontSize',20)
ylabel('R_{eco} (\mumol CO_2 m^-^2s^-^1)')
ax2.YAxisLocation = 'right';

z.TileSpacing = 'compact';
z.Padding = 'compact';
%% Compare concurrent measurements of PCO2 & GPP versus PCO2 & Reco

ind = find(ismember(datenum(USALQStreamDaily.Date_Time),datenum(USALQ.TIMESTAMP_END)));
USALQStreamDaily = USALQStreamDaily(ind,:);

ix = find(ismember(datenum(USALQ.TIMESTAMP_END),datenum(USALQStreamDaily.Date_Time)));
USALQ= USALQ(ix,:);

ind = find(ismember(datetime(year(USEDN.DateTime),month(USEDN.DateTime),day(USEDN.DateTime)),datetime(year(USEDNPCO2DAY.timestamp_sel),month(USEDNPCO2DAY.timestamp_sel),day(USEDNPCO2DAY.timestamp_sel))));
USEDN = USEDN(ind,:);

ix = find(ismember(datetime(year(USEDNPCO2DAY.timestamp_sel),month(USEDNPCO2DAY.timestamp_sel),day(USEDNPCO2DAY.timestamp_sel)),datetime(year(USEDN.DateTime),month(USEDN.DateTime),day(USEDN.DateTime))));
USEDNPCO2DAY= USEDNPCO2DAY(ix,:);

%% Statistics
% US-EDN SW
[R,P] = corrcoef(USEDNPCO2DAY.CO2(~isnan(USEDNPCO2DAY.CO2)&~isnan(USEDN.GPP_DT_uStar)),USEDN.GPP_DT_uStar(~isnan(USEDNPCO2DAY.CO2)&~isnan(USEDN.GPP_DT_uStar)))

[R,P] = corrcoef(USEDNPCO2DAY.CO2(~isnan(USEDNPCO2DAY.CO2)&~isnan(USEDN.Reco_DT_uStar)),USEDN.Reco_DT_uStar(~isnan(USEDNPCO2DAY.CO2)&~isnan(USEDN.Reco_DT_uStar)))

% USALQ PW
[R,P] = corrcoef(USALQ.pCO2(~isnan(USALQ.pCO2)&~isnan(USALQ.GPP_F)),USALQ.GPP_F(~isnan(USALQ.pCO2)&~isnan(USALQ.GPP_F)))

[R,P] = corrcoef(USALQ.pCO2(~isnan(USALQ.pCO2)&~isnan(USALQ.RECO_F)),USALQ.RECO_F(~isnan(USALQ.pCO2)&~isnan(USALQ.RECO_F)))

% USALQ SW
[R,P] = corrcoef(USALQStreamDaily.pCO2_ppm(~isnan(USALQStreamDaily.pCO2_ppm)&~isnan(USALQ.GPP_F)),USALQ.GPP_F(~isnan(USALQStreamDaily.pCO2_ppm)&~isnan(USALQ.GPP_F)))

[R,P] = corrcoef(USALQStreamDaily.pCO2_ppm(~isnan(USALQStreamDaily.pCO2_ppm)&~isnan(USALQ.Reco_F)),USALQ.Reco_F(~isnan(USALQStreamDaily.pCO2_ppm)&~isnan(USALQ.Reco_F)))

% SE-Deg PW
[R,P] = corrcoef(PorewaterSEDeg.PorewaterCO2(~isnan(PorewaterSEDeg.PorewaterCO2)&~isnan(PorewaterSEDeg.GPP_DT_uStar)),PorewaterSEDeg.GPP_DT_uStar(~isnan(PorewaterSEDeg.PorewaterCO2)&~isnan(PorewaterSEDeg.GPP_DT_uStar)))

[R,P] = corrcoef(PorewaterSEDeg.PorewaterCO2(~isnan(PorewaterSEDeg.PorewaterCO2)&~isnan(PorewaterSEDeg.Reco_DT_uStar)),PorewaterSEDeg.Reco_DT_uStar(~isnan(PorewaterSEDeg.PorewaterCO2)&~isnan(PorewaterSEDeg.Reco_DT_uStar)))

% SE-Deg SW
[R,P] = corrcoef(SEDeg.pCO2(~isnan(SEDeg.pCO2)&~isnan(SEDeg.GPP_DT_uStar)),SEDeg.GPP_DT_uStar(~isnan(SEDeg.pCO2)&~isnan(SEDeg.GPP_DT_uStar)))

[R,P] = corrcoef(SEDeg.pCO2(~isnan(SEDeg.pCO2)&~isnan(SEDeg.Reco_DT_uStar)),SEDeg.Reco_DT_uStar(~isnan(SEDeg.pCO2)&~isnan(SEDeg.Reco_DT_uStar)))
