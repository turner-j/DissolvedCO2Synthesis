% Clear workspace
close all;clear; clc;
run ImportUSEvMREddyProcOut.m
USGCE = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','US-GCE_DT');
USALQ = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','US-ALQ_DT');
run ImportALQstreampCO2.m
AllequashstreampCO2.Stage = 0.3048*AllequashstreampCO2.Stage; % convert stage from ft to m
USLos = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','US-Los_DT');
CADBB = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','CA-DBB');
SEDeg = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','SEDeg_DT');
LochVale = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','LochVale');
MBPPW1 = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','MBPPW1FCO2');
MBPPW2 = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','MBPPW2FCO2');

load('LochValepCO2.mat')
LochValeCO2 = Mastetal1998LochValeCO2flux.SaturatedAverageDailyCO2Fluxmmolm2d1;% Saturated soil FCO2
LochValeCO2 = LochValeCO2.*1000./86400;% Convert from units per day to seconds
load("LVDailyData.mat");LVDailyData = LVDailyData(LVDailyData.Site_ID==401723105400000,:);
run ImportLVMetData.m;LochValeClimate = WY1992to2019LochValeClimate((WY1992to2019LochValeClimate.station_name=='Andrews Creek weather station')&(WY1992to2019LochValeClimate.measurement_height==4),:);

addpath 'E:\SanDiskSecureAccess\all MATLAB files\Residence paper\USEvM'
addpath('E:\SanDiskSecureAccess\all MATLAB files\Residence paper\YKD')
YKDUB = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','YKDUB');
YKDB = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','YKDB');
run ImportBigCypresspCO2.m
USICs = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','USICs_DT');
USICs2 = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','USICs_DTinterp');
USOWC = readtable("E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls",'FileType','spreadsheet','Sheet','US-OWC_DT');
USOWC2 = readtable("E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls",'FileType','spreadsheet','Sheet','US-OWC_DT_interp');
USOWC = USOWC(6:end,:);
OWCsoilCO2 = USOWC(:,12:end);
USOWC.pCO2 = mean(table2array(OWCsoilCO2),2,'omitnan');% take the mean of each row

PLMdiel = readtable('E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls','Sheet','PLM');
PLM = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','PLM2');
PHMdiel = readtable('E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls','Sheet','PHM');
PHM = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','PHM2');

USHB1diel = readtable('E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls','Sheet','USHB1');
USHB1 = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','USHB1');

PorewaterSEDeg = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','PorewaterSEDeg_DT');
Troll = readtable("E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls",'FileType','spreadsheet','Sheet','Troll');
load USMyb.mat;USMyb2 = USMyb2(145:end,:);
APEX = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','APEX');
%% Calculate monthly average snow depth
ind = find(ismember(LochValeClimate.timestamp,LVDailyData.Date));
LochValeClimate = LochValeClimate(ind(1):ind(end),:);

s = string(LochValeClimate.Snow_depth);LochValeClimate.Snow_depth = double(s);
%% Eliminate unrealistic values
USLos.pCO2(USLos.pCO2<413)=NaN;
USALQ.pCO2(USALQ.pCO2<415)=NaN;
LochVale.pCO2(LochVale.pCO2<403.3)= NaN;
USMyb2.pCO2(USMyb2.pCO2<405.0)=NaN;
USMybREddyProcout.pCO2(USMybREddyProcout.pCO2<405.0)=NaN;
USMyb2 = USMyb2(145:end,:);

%% Alpine Sites
LVDiel = readtable("E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls",'FileType','spreadsheet','Sheet','LV_DT');
LVDiel.Snow_depth = str2double(LVDiel.Snow_depth);
ix = find(ismember(LVDiel.timestamp,LochValeClimate.timestamp));
LVDiel.Snow_depth = LochValeClimate.Snow_depth(ix,:);

LVGPP = [0;10;400;390;800;500;570;920;375;110];
LVGPP = (LVGPP.*1000)./86400;

DateStrings = {'1996-04-30';'1996-06-29';'1996-06-30';'1996-07-01';...
    '1996-07-10';'1996-07-11';'1996-07-25';'1996-08-01';'1996-09-01';'1996-10-20'};
dt = datetime(DateStrings,'InputFormat','yyyy-MM-dd');

monthlyFCO2 = accumarray(month(dt),LVGPP,[],@(x)mean(x,'omitnan'),NaN);
[~,id] = unique(month(dt));
gramC = sum(monthlyFCO2(~isnan(monthlyFCO2)).*(1E-6).*(12).*(2.628e+6));% yearly total GPP
disp(gramC)
SEM = std((monthlyFCO2(~isnan(monthlyFCO2)).*(1E-6).*(12).*(2.628e+6)))./sqrt(length(monthlyFCO2));
disp(SEM)

t = tiledlayout(2,2);
nexttile(1)
title(t,'Seasonal GPP Cycling','FontSize',20)
ylabel(t,'GPP (\mumol CO_2 m^-^2s^-^1)','FontSize',20)
xlabel(t,'Day of Year','FontSize',20)

title('Alpine')
hold on
err = accumarray(month(dt),LVGPP,[],@(x)std(x,'omitnan'),NaN);
errorbar(day(dt(id),'dayofyear'),monthlyFCO2(~isnan(monthlyFCO2)),err(~isnan(err)),'o-','Color',[0.2118 0.3529 0.4784],'MarkerFaceColor',[0.2118 0.3529 0.4784],'MarkerEdgeColor',[0.2118 0.3529 0.4784],'LineWidth',1.5,'MarkerSize',6)
hold off
xlim([0 370])
xticks(0:50:350)
xticklabels({''})
set(gca,'FontSize',20)
box on
xline(214);

% Prairie Potholes and Karsts
nexttile(3)
title('Prairie Potholes and Karsts')
hold on
plot(1,nan,'o-','Color',[0.2118 0.3529 0.4784],'MarkerFaceColor',[0.2118 0.3529 0.4784],'MarkerEdgeColor',[0.2118 0.3529 0.4784],'LineWidth',1.5,'MarkerSize',6)

monthlyFCO2 = accumarray(month(MBPPW1.DATE),MBPPW1.GPP_f,[],@(x)mean(x,'omitnan'),NaN);
[~,id] = unique(month(MBPPW1.DATE));
err = accumarray(month(MBPPW1.DATE),MBPPW1.GPP_f,[],@(x)std(x,'omitnan'),NaN);
errorbar(day(MBPPW1.DATE(id),'dayofyear'),monthlyFCO2(~isnan(monthlyFCO2)),err(~isnan(err)),'o-','Color',[0.0588 0.4902 0.5216],'MarkerFaceColor',[0.0588 0.4902 0.5216],'MarkerEdgeColor',[0.0588 0.4902 0.5216],'LineWidth',1.5,'MarkerSize',6)

gramC = sum(monthlyFCO2(~isnan(monthlyFCO2)).*(1E-6).*(12).*(2.628e+6));% yearly total GPP
disp(gramC)
SEM = std((monthlyFCO2(~isnan(monthlyFCO2)).*(1E-6).*(12).*(2.628e+6)))./sqrt(length(monthlyFCO2));
disp(SEM)

monthlyFCO2 = accumarray(month(MBPPW2.DATE),MBPPW2.GPP_f,[],@(x)mean(x,'omitnan'),NaN);
[~,id] = unique(month(MBPPW2.DATE));
err = accumarray(month(MBPPW2.DATE),MBPPW2.GPP_f,[],@(x)std(x,'omitnan'),NaN);
errorbar(day(MBPPW2.DATE(id),'dayofyear'),monthlyFCO2(~isnan(monthlyFCO2)),err(~isnan(err)),'o-','Color',[0.5961 0.8706 0.8902],'MarkerFaceColor',[0.5961 0.8706 0.8902],'MarkerEdgeColor',[0.5961 0.8706 0.8902],'LineWidth',1.5,'MarkerSize',6)

gramC = sum(monthlyFCO2(~isnan(monthlyFCO2)).*(1E-6).*(12).*(2.628e+6));% yearly total GPP
disp(gramC)
SEM = std((monthlyFCO2(~isnan(monthlyFCO2)).*(1E-6).*(12).*(2.628e+6)))./sqrt(length(monthlyFCO2));
disp(SEM)

plot(1,nan,'o-','Color',[0.0588 0.3216 0.1294],'MarkerFaceColor',[0.0588 0.3216 0.1294],'MarkerEdgeColor',[0.0588 0.3216 0.1294],'LineWidth',1.5,'MarkerSize',6)
plot(1,nan,'o--k','MarkerFaceColor',[0.2039 0.6196 0.3137],'MarkerEdgeColor','k','LineWidth',1.5,'MarkerSize',5)
plot(1,nan,'o-','Color',[0.3686 1.0000 0.4941],'MarkerFaceColor',[0.3686 1.0000 0.4941],'MarkerEdgeColor',[0.3686 1.0000 0.4941],'LineWidth',1.5,'MarkerSize',6)
plot(1,nan,'o-','Color',[0.6706 0.7882 0.7176],'MarkerFaceColor',[0.6706 0.7882 0.7176],'MarkerEdgeColor',[0.6706 0.7882 0.7176],'LineWidth',1.5,'MarkerSize',6)
plot(1,nan,'o-','Color',[0.6627 0.9882 0.0118],'MarkerFaceColor',[0.6627 0.9882 0.0118],'MarkerEdgeColor','k','LineWidth',1.5,'MarkerSize',5)
plot(1,nan,'o-','Color',[0.3333 0.4902 0.0235],'MarkerFaceColor',[0.3333 0.4902 0.0235],'MarkerEdgeColor',[0.3333 0.4902 0.0235],'LineWidth',1.5,'MarkerSize',6)
plot(1,nan,'o-k','MarkerFaceColor',[0.7373 0.8196 0.2118],'MarkerEdgeColor',[0.7373 0.8196 0.2118],'LineWidth',1.5,'MarkerSize',5)

hold off
xlim([0 370])
xline(186);
set(gca,'FontSize',20)
box on
leg = legend({'LV','MBPPW1','MBPPW2','SE-Deg','CA-DBB','Troll','US-Myb','US-Los (PW)','US-ALQ (PW)','SE-Deg (PW)'},...
    'FontSize',20,'Orientation','horizontal');
leg.Layout.Tile = 'south';

% fens and bogs
run ImportTrollFCO2.m
run ImportTrollpCO2.m

% Plot pCO2
ax2 = nexttile(2);
hold on
monthlyFCO2 = accumarray(month(SEDeg.DateTime),SEDeg.GPP_DT_uStar,[],@(x)mean(x,'omitnan'),NaN);
[~,id] = unique(month(SEDeg.DateTime));
err = accumarray(month(SEDeg.DateTime),SEDeg.GPP_DT_uStar,[],@(x)std(x,'omitnan'),NaN);
errorbar(day(SEDeg.DateTime(id),'dayofyear'),monthlyFCO2(~isnan(monthlyFCO2)),err(~isnan(err)),'o-','Color',[0.0588 0.3216 0.1294],'MarkerFaceColor',[0.0588 0.3216 0.1294],'MarkerEdgeColor',[0.0588 0.3216 0.1294],'LineWidth',1.5,'MarkerSize',6)

gramC = sum(monthlyFCO2.*(1E-6).*(12).*(2.628e+6));
disp(gramC)% yearly total GPP
SEM = std((monthlyFCO2(~isnan(monthlyFCO2)).*(1E-6).*(12).*(2.628e+6)))./sqrt(length(monthlyFCO2));
disp(SEM)

monthlyFCO2 = accumarray(month(CADBB.DATE),CADBB.GPP,[],@(x)mean(x,'omitnan'),NaN);
[~,id] = unique(month(CADBB.DATE));
err = accumarray(month(CADBB.DATE),CADBB.GPP,[],@(x)std(x,'omitnan'),NaN);
errorbar(day(CADBB.DATE(id),'dayofyear'),monthlyFCO2(~isnan(monthlyFCO2)),err(~isnan(err)),'o--k','MarkerFaceColor',[0.2039 0.6196 0.3137],'MarkerEdgeColor','k','LineWidth',1.5,'MarkerSize',6)

gramC = sum(monthlyFCO2.*(1E-6).*(12).*(2.628e+6));% yearly total GPP
disp(gramC)
SEM = std((monthlyFCO2(~isnan(monthlyFCO2)).*(1E-6).*(12).*(2.628e+6)))./sqrt(length(monthlyFCO2));
disp(SEM)

monthlyFCO2 = accumarray(month(TrollFCO2.DateTime),TrollFCO2.GPP_uStar_f,[],@(x)mean(x,'omitnan'),NaN);
[~,id] = unique(month(TrollFCO2.DateTime));
err = accumarray(month(TrollFCO2.DateTime),TrollFCO2.GPP_uStar_f,[],@(x)std(x,'omitnan'),NaN);
errorbar(day(TrollFCO2.DateTime(id),'dayofyear'),monthlyFCO2(~isnan(monthlyFCO2)),err(~isnan(err)),'o-','Color',[0.3686 1.0000 0.4941],'MarkerFaceColor',[0.3686 1.0000 0.4941],'MarkerEdgeColor',[0.3686 1.0000 0.4941],'LineWidth',1.5,'MarkerSize',6)

gramC = sum(monthlyFCO2.*(1E-6).*(12).*(2.628e+6));% yearly total GPP
disp(gramC)
SEM = std((monthlyFCO2(~isnan(monthlyFCO2)).*(1E-6).*(12).*(2.628e+6)))./sqrt(length(monthlyFCO2));
disp(SEM)

monthlyFCO2 = accumarray(month(USMyb2.DateTime),USMyb2.GPP_DT_uStar,[],@(x)mean(x,'omitnan'),NaN);
[~,id] = unique(month(USMyb2.DateTime));
err = accumarray(month(USMyb2.DateTime),USMyb2.GPP_DT_uStar,[],@(x)std(x,'omitnan'),NaN);
errorbar(day(USMyb2.DateTime(id),'dayofyear'),monthlyFCO2,err,'o-','Color',[0.6706 0.7882 0.7176],'MarkerFaceColor',[0.6706 0.7882 0.7176],'MarkerEdgeColor',[0.6706 0.7882 0.7176],'LineWidth',1.5,'MarkerSize',6)

gramC = sum(monthlyFCO2.*(1E-6).*(12).*(2.628e+6));
disp(gramC)% yearly total GPP
SEM = std((monthlyFCO2(~isnan(monthlyFCO2)).*(1E-6).*(12).*(2.628e+6)))./sqrt(length(monthlyFCO2));
disp(SEM)

xline(181);
xticklabels({''})
xticks(0:50:350)
hold off
xlim([0 370])
box on
title('Fens, Bogs, and Marshes')
set(gca,'FontSize',20)
ax2.YAxisLocation = 'right';

% Porewater
load SEDegPorewater.mat
USLos2 = readtable("E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls",'FileType','spreadsheet','Sheet','US-Los');
USALQ2 = readtable("E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls",'FileType','spreadsheet','Sheet','US-ALQ_DT');
SEDegREddyProcOut2 = SEDegREddyProcOut2(7915:77256,:);

load('USLosFCO2_2020.mat')

ax = nexttile(4);
title('Porewater')
hold on
monthlyFCO2 = accumarray(month(NEELos.TIMESTAMP_END),NEELos.GPP_PI_F,[],@(x)mean(x,'omitnan'),NaN);
[~,id] = unique(month(NEELos.TIMESTAMP_END));
err = accumarray(month(NEELos.TIMESTAMP_END),NEELos.GPP_PI_F,[],@(x)std(x,'omitnan'),NaN);
errorbar(day(NEELos.TIMESTAMP_END(id),'dayofyear'),monthlyFCO2(~isnan(monthlyFCO2)),err(~isnan(err)),'o-','Color',[0.6627 0.9882 0.0118],'MarkerFaceColor',[0.6627 0.9882 0.0118],'MarkerEdgeColor','k','LineWidth',1.5,'MarkerSize',6)

gramC = sum(monthlyFCO2.*(1E-6).*(12).*(2.628e+6));
disp(gramC)% yearly total GPP
SEM = std((monthlyFCO2(~isnan(monthlyFCO2)).*(1E-6).*(12).*(2.628e+6)))./sqrt(length(monthlyFCO2));
disp(SEM)

monthlyFCO2 = accumarray(month(USALQ.TIMESTAMP_END),USALQ.GPP_F,[],@(x)mean(x,'omitnan'),NaN);
[~,id] = unique(month(USALQ.TIMESTAMP_END));
err = accumarray(month(USALQ.TIMESTAMP_END),USALQ.GPP_F,[],@(x)std(x,'omitnan'),NaN);
errorbar(day(USALQ.TIMESTAMP_END(id),'dayofyear'),monthlyFCO2(~isnan(monthlyFCO2)),err(~isnan(err)),'o-','Color',[0.3333 0.4902 0.0235],'MarkerFaceColor',[0.3333 0.4902 0.0235],'MarkerEdgeColor',[0.3333 0.4902 0.0235],'LineWidth',1.5,'MarkerSize',6)

gramC = sum(monthlyFCO2.*(1E-6).*(12).*(2.628e+6));
disp(gramC)% yearly total GPP
SEM = std((monthlyFCO2(~isnan(monthlyFCO2)).*(1E-6).*(12).*(2.628e+6)))./sqrt(length(monthlyFCO2));
disp(SEM)

monthlyFCO2 = accumarray(month(Deg3.DateTime),Deg3.GPP_DT_uStar,[],@(x)mean(x,'omitnan'),NaN);
[~,id] = unique(month(Deg3.DateTime));
err = accumarray(month(Deg3.DateTime),Deg3.GPP_DT_uStar,[],@(x)std(x,'omitnan'),NaN);
errorbar(day(Deg3.DateTime(id),'dayofyear'),monthlyFCO2(~isnan(monthlyFCO2)),err(~isnan(err)),'o-k','MarkerFaceColor',[0.7373 0.8196 0.2118],'MarkerEdgeColor',[0.7373 0.8196 0.2118],'LineWidth',1.5,'MarkerSize',6)

box on
ax.YAxisLocation = 'right';
t.TileSpacing = 'compact';
t.Padding = 'compact';
xline(182);
xlim([0 370])
set(gca,'FontSize',20)
%% Seasonal cycles of pCO2
t = tiledlayout(2,2);
nexttile(1)
title(t,'Seasonal P_{CO2} Cycling','FontSize',20)
ylabel(t,'CO_2({\itaq}) (ppmv)','FontSize',20)
xlabel(t,'Day of Year','FontSize',20)

title('Alpine')
hold on
monthlypCO2 = accumarray(month(LVDailyData.Date),LVDailyData.Dissolved_CO2_ppmv,[],@(x)mean(x,'omitnan'),NaN);
[~,id] = unique(month(LVDailyData.Date));
err = accumarray(month(LVDailyData.Date),LVDailyData.Dissolved_CO2_ppmv,[],@(x)std(x,'omitnan'),NaN);
errorbar(day(LVDailyData.Date(id),'dayofyear'),monthlypCO2,err,'o-','Color',[0.2118 0.3529 0.4784],'MarkerFaceColor',[0.2118 0.3529 0.4784],'MarkerEdgeColor',[0.2118 0.3529 0.4784],'LineWidth',1.5,'MarkerSize',6)
hold off
xlim([0 370])
xticks(0:50:350)
xticklabels({''})
set(gca,'FontSize',20)
box on
xline(92);

% Prairie Potholes and Karsts
nexttile(3)
title('Prairie Potholes and Karsts')
hold on
plot(1,nan,'o-','Color',[0.2118 0.3529 0.4784],'MarkerFaceColor',[0.2118 0.3529 0.4784],'MarkerEdgeColor',[0.2118 0.3529 0.4784],'LineWidth',1.5,'MarkerSize',6)

addpath('E:\SanDiskSecureAccess\all MATLAB files\Residence paper\PrairiePotholes')
run ImportDissolvedGas.m
% Sort data by site
SITEID = categories(dissolvedgas.WETLAND_ID);
myTables = [];
fn = regexprep(SITEID, '\s', '');% get rid of spaces in names
fn = regexprep(fn, '\W', '');

for i = 1:length(SITEID)
    tablenow = dissolvedgas((dissolvedgas.WETLAND_ID==SITEID(i)),:);
    tablenow.DATE = tablenow.DATE + timeofday(tablenow.TIME); 
    yr = year(tablenow.DATE);
    idx = find(yr~=2014);
    pCO2Tables.(fn{i}) = tablenow(idx,:);
end
clearvars SITEID tablenow

count = 0;

% Delete flux tables without pCO2 data and create daily averages 
% Plot monthly averages 
for n = 1:length(fn)
    Tableout = pCO2Tables.(fn{n});
    mo = unique(month(Tableout.DATE));
    [~,~,BC] = groupsummary(Tableout.CO2_CONC,(month(Tableout.DATE)),'mean');
    if (sum(BC>=10)>=5)% only display sites with more than 5 measurements per month for ten months
    monthlypCO2 = accumarray(month(Tableout.DATE),Tableout.CO2_CONC,[],@(x)mean(x,'omitnan'),NaN);
    [~,id] = unique(month(Tableout.DATE));
    err = accumarray(month(Tableout.DATE),Tableout.CO2_CONC,[],@(x)std(x,'omitnan'),NaN);
    errorbar(day(Tableout.DATE(id),'dayofyear'),monthlypCO2(~isnan(monthlypCO2)),err(~isnan(err)),'o-','Color',[0.5882 0.7569 0.8314],'MarkerFaceColor',[0.7569 0.8980 0.9608],'MarkerEdgeColor',[0.7569 0.8980 0.9608],'LineWidth',1.5,'MarkerSize',6)
    count = count+1;
    range(monthlypCO2)
    end
end

monthlypCO2 = accumarray(month(MBPPW1.DATE),MBPPW1.pCO2,[],@(x)mean(x,'omitnan'),NaN);
[~,id] = unique(month(MBPPW1.DATE));
err = accumarray(month(MBPPW1.DATE),MBPPW1.pCO2,[],@(x)std(x,'omitnan'),NaN);
errorbar(day(MBPPW1.DATE(id),'dayofyear'),monthlypCO2(~isnan(monthlypCO2)),err(~isnan(err)),'o-','Color',[0.0588 0.4902 0.5216],'MarkerFaceColor',[0.0588 0.4902 0.5216],'MarkerEdgeColor',[0.0588 0.4902 0.5216],'LineWidth',1.5,'MarkerSize',6)

monthlypCO2 = accumarray(month(MBPPW2.DATE),MBPPW2.pCO2,[],@(x)mean(x,'omitnan'),NaN);
[~,id] = unique(month(MBPPW2.DATE));
err = accumarray(month(MBPPW2.DATE),MBPPW2.pCO2,[],@(x)std(x,'omitnan'),NaN);
errorbar(day(MBPPW2.DATE(id),'dayofyear'),monthlypCO2(~isnan(monthlypCO2)),err(~isnan(err)),'o-k','MarkerFaceColor',[0.5961 0.8706 0.8902],'MarkerEdgeColor','k','LineWidth',1.5,'MarkerSize',6)

plot(1,nan,'o-','Color',[0.0588 0.3216 0.1294],'MarkerFaceColor',[0.0588 0.3216 0.1294],'MarkerEdgeColor',[0.0588 0.3216 0.1294],'LineWidth',1.5,'MarkerSize',6)
plot(1,nan,'o--k','MarkerFaceColor',[0.2039 0.6196 0.3137],'MarkerEdgeColor','k','LineWidth',1.5,'MarkerSize',5)
plot(1,nan,'o-','Color',[0.3686 1.0000 0.4941],'MarkerFaceColor',[0.3686 1.0000 0.4941],'MarkerEdgeColor',[0.3686 1.0000 0.4941],'LineWidth',1.5,'MarkerSize',6)
plot(1,nan,'o-','Color',[0.6706 0.7882 0.7176],'MarkerFaceColor',[0.6706 0.7882 0.7176],'MarkerEdgeColor',[0.6706 0.7882 0.7176],'LineWidth',1.5,'MarkerSize',6)
plot(1,nan,'o-','Color',[0.6627 0.9882 0.0118],'MarkerFaceColor',[0.6627 0.9882 0.0118],'MarkerEdgeColor','k','LineWidth',1.5,'MarkerSize',5)
plot(1,nan,'o-','Color',[0.3333 0.4902 0.0235],'MarkerFaceColor',[0.3333 0.4902 0.0235],'MarkerEdgeColor',[0.3333 0.4902 0.0235],'LineWidth',1.5,'MarkerSize',6)
plot(1,nan,'o-','Color',[0.1451 0.2118 0.0078],'MarkerFaceColor',[0.1451 0.2118 0.0078],'MarkerEdgeColor',[0.1451 0.2118 0.0078],'LineWidth',1.5,'MarkerSize',6)
plot(1,nan,'o-k','MarkerFaceColor',[0.7373 0.8196 0.2118],'MarkerEdgeColor',[0.7373 0.8196 0.2118],'LineWidth',1.5,'MarkerSize',5)

hold off
xlim([0 370])
xline(209);

set(gca,'FontSize',20)
box on
leg = legend({'LV','Cotton.','','','','','','','MBPPW1','MBPPW2','SE-Deg','CA-DBB','Troll.','US-Myb','US-Los (PW)','US-ALQ (PW)','US-OWC (PW)','SE-Deg (PW)'},...
    'FontSize',20,'Orientation','horizontal','NumColumns',6);
leg.Layout.Tile = 'south';

% fens and bogs
% Plot pCO2
ax2 = nexttile(2);
hold on
monthlypCO2 = accumarray(month(SEDeg.DateTime),SEDeg.pCO2,[],@(x)mean(x,'omitnan'),NaN);
[~,id] = unique(month(SEDeg.DateTime));
err = accumarray(month(SEDeg.DateTime),SEDeg.pCO2,[],@(x)std(x,'omitnan'),NaN);
errorbar(day(SEDeg.DateTime(id),'dayofyear'),monthlypCO2(~isnan(monthlypCO2)),err(~isnan(err)),'o-','Color',[0.0588 0.3216 0.1294],'MarkerFaceColor',[0.0588 0.3216 0.1294],'MarkerEdgeColor',[0.0588 0.3216 0.1294],'LineWidth',1.5,'MarkerSize',6)

monthlypCO2 = accumarray(month(CADBB.DATE),CADBB.pCO2,[],@(x)mean(x,'omitnan'),NaN);
[~,id] = unique(month(CADBB.DATE));
id = id(~isnan(monthlypCO2));
err = accumarray(month(CADBB.DATE),CADBB.pCO2,[],@(x)std(x,'omitnan'),NaN);
errorbar(day(CADBB.DATE(id),'dayofyear'),monthlypCO2(~isnan(monthlypCO2)),err(~isnan(err)),'o--k','MarkerFaceColor',[0.2039 0.6196 0.3137],'MarkerEdgeColor','k','LineWidth',1.5,'MarkerSize',6)

monthlypCO2 = accumarray(month(TrollpCO2.FieldDate),TrollpCO2.pCO2atm,[],@(x)mean(x,'omitnan'),NaN);
range(monthlypCO2)
[~,id] = unique(month(TrollpCO2.FieldDate));
err = accumarray(month(TrollpCO2.FieldDate),TrollpCO2.pCO2atm,[],@(x)std(x,'omitnan'),NaN);
errorbar(day(TrollpCO2.FieldDate(id(~isnan(monthlypCO2))),'dayofyear'),monthlypCO2(~isnan(monthlypCO2)),err(~isnan(err)),'o-','Color',[0.3686 1.0000 0.4941],'MarkerFaceColor',[0.3686 1.0000 0.4941],'MarkerEdgeColor',[0.3686 1.0000 0.4941],'LineWidth',1.5,'MarkerSize',6)

monthlypCO2 = accumarray(month(USMyb2.DateTime),USMyb2.pCO2,[],@(x)mean(x,'omitnan'),NaN);
[~,id] = unique(month(USMyb2.DateTime));
err = accumarray(month(USMyb2.DateTime),USMyb2.pCO2,[],@(x)std(x,'omitnan'),NaN);
errorbar(day(USMyb2.DateTime(id),'dayofyear'),monthlypCO2,err,'o-','Color',[0.6706 0.7882 0.7176],'MarkerFaceColor',[0.6706 0.7882 0.7176],'MarkerEdgeColor',[0.6706 0.7882 0.7176],'LineWidth',1.5,'MarkerSize',6)

monthlypCO2 = accumarray(month(AllequashstreampCO2.VarName1),AllequashstreampCO2.CO5,[],@(x)mean(x,'omitnan'),NaN);
[~,id] = unique(month(AllequashstreampCO2.VarName1));
err = accumarray(month(AllequashstreampCO2.VarName1),AllequashstreampCO2.CO5,[],@(x)std(x,'omitnan'),NaN);
errorbar(day(AllequashstreampCO2.VarName1(id),'dayofyear'),monthlypCO2(~isnan(monthlypCO2)),err(~isnan(err)),'o-k','MarkerFaceColor',[0.5608 0.6588 0.5843],'MarkerEdgeColor','k','LineWidth',1.5,'MarkerSize',6)
range(monthlypCO2)

xline(222.4);
xticklabels({''})
xticks(0:50:350)
hold off
xlim([0 370])
box on
title('Fens, Bogs, and Marshes')
set(gca,'FontSize',20)
ax2.YAxisLocation = 'right';

% Porewater
ax = nexttile(4);
title('Porewater')
hold on
monthlypCO2 = accumarray(month(USLos.TIMESTAMP_END),USLos.pCO2,[],@(x)mean(x,'omitnan'),NaN);
[~,id] = unique(month(USLos.TIMESTAMP_END));
err = accumarray(month(USLos.TIMESTAMP_END),USLos.pCO2,[],@(x)std(x,'omitnan'),NaN);
errorbar(day(USLos.TIMESTAMP_END(id),'dayofyear'),monthlypCO2(~isnan(monthlypCO2)),err(~isnan(err)),'o-','Color',[0.6627 0.9882 0.0118],'MarkerFaceColor',[0.6627 0.9882 0.0118],'MarkerEdgeColor','k','LineWidth',1.5,'MarkerSize',6)

monthlypCO2 = accumarray(month(USALQ.TIMESTAMP_END),USALQ.pCO2,[],@(x)mean(x,'omitnan'),0);
[~,id] = unique(month(USALQ.TIMESTAMP_END));
err = accumarray(month(USALQ.TIMESTAMP_END),USALQ.pCO2,[],@(x)std(x,'omitnan'),NaN);
errorbar(day(USALQ.TIMESTAMP_END(id),'dayofyear'),monthlypCO2(monthlypCO2~=0),err(monthlypCO2~=0),'o-','Color',[0.3333 0.4902 0.0235],'MarkerFaceColor',[0.3333 0.4902 0.0235],'MarkerEdgeColor',[0.3333 0.4902 0.0235],'LineWidth',1.5,'MarkerSize',6)

monthlypCO2 = accumarray(month(USOWC.TIMESTAMP_END),USOWC.pCO2,[],@(x)mean(x,'omitnan'),NaN);
[~,id] = unique(month(USOWC.TIMESTAMP_END));
err = accumarray(month(USOWC.TIMESTAMP_END),USOWC.pCO2,[],@(x)std(x,'omitnan'),NaN);
errorbar(day(USOWC.TIMESTAMP_END(id),'dayofyear'),monthlypCO2(~isnan(monthlypCO2)),err(~isnan(err)),'o-','Color',[0.1451 0.2118 0.0078],'MarkerFaceColor',[0.1451 0.2118 0.0078],'MarkerEdgeColor',[0.1451 0.2118 0.0078],'LineWidth',1.5,'MarkerSize',6)

monthlypCO2 = accumarray(month(Deg3.DateTime),Deg3.PorewaterCO2,[],@(x)mean(x,'omitnan'),NaN);
[~,id] = unique(month(Deg3.DateTime));
err = accumarray(month(Deg3.DateTime),Deg3.PorewaterCO2,[],@(x)std(x,'omitnan'),NaN);
errorbar(day(Deg3.DateTime(id),'dayofyear'),monthlypCO2(~isnan(monthlypCO2)),err(~isnan(err)),'o-k','MarkerFaceColor',[0.7373 0.8196 0.2118],'MarkerEdgeColor',[0.7373 0.8196 0.2118],'LineWidth',1.5,'MarkerSize',6)

box on
ax.YAxisLocation = 'right';
t.TileSpacing = 'compact';
t.Padding = 'compact';
xline(223.5);
xlim([0 370])
set(gca,'FontSize',20)

