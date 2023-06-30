% Clear workspace
close all;clear; clc;
run ImportUSEvMREddyProcOut.m
USGCE = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','US-GCE_DT');
USALQ = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','US-ALQ_DT');
load USALQStreamDaily.mat
ind = find(ismember(USALQStreamDaily.Date_Time,USALQ.TIMESTAMP_END));
USALQStreamDaily = USALQStreamDaily(ind,:);

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

load SEDegPorewater.mat
USEDN = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','US-EDN_DT');
load SouthPennines.mat
%% Calculate monthly average snow depth
ind = find(ismember(LochValeClimate.timestamp,LVDailyData.Date));
LochValeClimate = LochValeClimate(ind(1):ind(end),:);

s = string(LochValeClimate.Snow_depth);LochValeClimate.Snow_depth = double(s);
monthlysnow = accumarray(month(LochValeClimate.timestamp),LochValeClimate.Snow_depth,[],@(x)mean(x,'omitnan'),NaN);
%% Eliminate unrealistic values
USLos.pCO2(USLos.pCO2<413)=NaN;
USALQ.pCO2(USALQ.pCO2<415)=NaN;
LochVale.pCO2(LochVale.pCO2<403.3)= NaN;
USMyb2.pCO2(USMyb2.pCO2<405.0)=NaN;
USMybREddyProcout.pCO2(USMybREddyProcout.pCO2<405.0)=NaN;

LVDiel = readtable("E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls",'FileType','spreadsheet','Sheet','LV_DT');
LVDiel.Snow_depth = str2double(LVDiel.Snow_depth);
ix = find(ismember(LVDiel.timestamp,LochValeClimate.timestamp));
LVDiel.Snow_depth = LochValeClimate.Snow_depth(ix,:);

% fens and bogs
run ImportTrollFCO2.m
run ImportTrollpCO2.m

%% bin average according to air temperature
binEdges = 0:2:30;

groupTA3 = discretize(SEDeg.Tair_f,binEdges);
[meanTA3,BG3] = groupsummary(SEDeg.pCO2,groupTA3,'mean');
[R,P] = corrcoef(binEdges(BG3(~isnan(BG3))),meanTA3(~isnan(BG3)));
disp([R,P])

groupTAC = discretize(Deg3.Tair_f,binEdges);
[meanTAC,BGC] = groupsummary(Deg3.PorewaterCO2,groupTAC,'mean');
[R,P] = corrcoef(binEdges(BGC(~isnan(BGC))),meanTAC(~isnan(BGC)));
disp([R,P])

groupTA4 = discretize(USLos.Tair_f,binEdges);
[meanTA4,BG4] = groupsummary(USLos.pCO2,groupTA4,'mean');
[R,P] = corrcoef(binEdges(BG4(~isnan(BG4))),meanTA4(~isnan(BG4)));
disp([R,P])

groupTA = discretize(USALQ.TA_F,binEdges);
[meanTA,BG] = groupsummary(USALQ.pCO2,groupTA,'mean');
[R,P] = corrcoef(binEdges(BG(~isnan(meanTA)&~isnan(BG))),meanTA(~isnan(meanTA)&~isnan(BG)));
disp([R,P])

groupTA2 = discretize(CADBB.Air_temp_2m,binEdges);
[meanTA2,BG2] = groupsummary(CADBB.pCO2,groupTA2,'mean');
[R,P] = corrcoef(binEdges(BG2(~isnan(meanTA2)&~isnan(BG2))),meanTA2(~isnan(meanTA2)&~isnan(BG2)));
disp([R,P])

groupTA6 = discretize(LochVale.T_air,binEdges);
[meanTA6,BG6] = groupsummary(LochVale.pCO2,groupTA6,'mean');
[R,P] = corrcoef(binEdges(BG6(~isnan(BG6))),meanTA6(~isnan(BG6)));
disp([R,P])

groupTA9 = discretize(MBPPW1.TA_1_1_1,binEdges);
[meanTA9,BG9] = groupsummary(MBPPW1.pCO2,groupTA9,'mean');
[R,P] = corrcoef(binEdges(BG9(~isnan(BG9))),meanTA9(~isnan(BG9)));
disp([R,P])

groupTA1 = discretize(USICs.Tair_f,binEdges);
[meanTA1,BG1] = groupsummary(USICs.pCO2,groupTA1,'mean');
[R,P] = corrcoef(binEdges(BG1(~isnan(BG1))),meanTA1(~isnan(BG1)));
disp([R,P])

groupTAB = discretize(USMyb2.Tair_f,binEdges);
[meanTAB,BGB] = groupsummary(USMyb2.pCO2,groupTAB,'mean');
[R,P] = corrcoef(binEdges(BGB(~isnan(BGB)&~isnan(meanTAB))),meanTAB(~isnan(BGB)&~isnan(meanTAB)));
disp([R,P])

groupTA7 = discretize(Troll.Tair,binEdges);
[meanTA7,BG7] = groupsummary(Troll.pCO2,groupTA7,'mean');
[R,P] = corrcoef(binEdges(BG7(~isnan(BG7)&~isnan(meanTA7))),meanTA7(~isnan(BG7)&~isnan(meanTA7)));
disp([R,P])

groupTA8 = discretize(USOWC.TA,binEdges);
[meanTA8,BG8] = groupsummary(USOWC.pCO2,groupTA8,'mean');
[R,P] = corrcoef(binEdges(BG8(~isnan(BG8))),meanTA8(~isnan(BG8)));
disp([R,P])

groupTAM = discretize(PHM.TA_F,binEdges);
[meanTAM,BGM] = groupsummary(PHM.pCO2,groupTAM,'mean');
[R,P] = corrcoef(binEdges(BGM(~isnan(BGM))),meanTAM(~isnan(BGM)));
disp([R,P])

groupTAX = discretize(USEDN.Tair_f,binEdges);
[meanTAX,BGX] = groupsummary(USEDN.CO2,groupTAX,'mean');
[R,P] = corrcoef(binEdges(BGX(~isnan(BGX))),meanTAX(~isnan(BGX)));
disp([R,P])

% now for the cotton. sites

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

% Prairie Potholes

for i = 1:length(fn)
    tablenow3 = pCO2Tables.(fn{(i)});
    %Take the averages across each day for all variables
    for k= 4:width(tablenow3)
        tablenow3.(k)(tablenow3.(k)==-9999) = NaN;       
    end
    newpCO2Tables.(fn{i}) = tablenow3;
end

for i = 1:length(fn)
    Tableout = newpCO2Tables.(fn{i});
        groupTA11 = discretize(Tableout.AIR_TEMP_C,binEdges);
        [meanTA11,BGZ] = groupsummary(Tableout.CO2_CONC,groupTA11,'mean');
        if length(BGZ(~isnan(BGZ)))>=3
            [R,P] = corrcoef(binEdges(BGZ(~isnan(BGZ))),meanTA11(~isnan(BGZ)));
            disp([R,P])
        end
end

%% Plotting

t = tiledlayout(3,2);
nexttile(1)
title(t,'P_{CO2} vs. T_{air}','FontSize',20)
plot(binEdges(BG6(~isnan(BG6))),meanTA6(~isnan(BG6)),'o--','Color',[0.2118 0.3529 0.4784],'LineWidth',2,'MarkerFaceColor',[0.2118 0.3529 0.4784],'MarkerEdgeColor',[0.2118 0.3529 0.4784])% LV
title('Alpine')
xlim([0 30])
box on
ylim([0 5000])
set(gca,'FontSize',20)
xticklabels('')

ax = nexttile(2);
title('Fens, Bogs, and Marshes')
hold on
plot(binEdges(BG3(~isnan(BG3))),meanTA3(~isnan(BG3)),'o-','LineWidth',2,'Color',[0.0588 0.3216 0.1294],'MarkerFaceColor',[0.0588 0.3216 0.1294],'MarkerEdgeColor',[0.0588 0.3216 0.1294])% SE-Deg
plot(binEdges(BG2(~isnan(meanTA2)&~isnan(BG2))),meanTA2(~isnan(meanTA2)&~isnan(BG2)),'o--','LineWidth',2,'MarkerFaceColor',[0.2039 0.6196 0.3137],'MarkerEdgeColor',[0.2039 0.6196 0.3137],'Color',[0.2039 0.6196 0.3137])% CA-DBB
plot(binEdges(BG1(~isnan(BG1))),meanTA1(~isnan(BG1)),'o--','LineWidth',2,'Color',[0.3647 0.7608 0.5765],'MarkerFaceColor',[0.3647 0.7608 0.5765],'MarkerEdgeColor',[0.3647 0.7608 0.5765])% US-ICs
plot(binEdges(BGB(~isnan(BGB)&~isnan(meanTAB))),meanTAB(~isnan(BGB)&~isnan(meanTAB)),'o--','LineWidth',2,'Color',[0.6706 0.7882 0.7176],'MarkerFaceColor',[0.6706 0.7882 0.7176],'MarkerEdgeColor',[0.6706 0.7882 0.7176])% US-Myb
plot(binEdges(BG7(~isnan(BG7)&~isnan(meanTA7))),meanTA7(~isnan(BG7)&~isnan(meanTA7)),'o--','LineWidth',2,'MarkerFaceColor',[0.3686 1.0000 0.4941],'MarkerEdgeColor',[0.3686 1.0000 0.4941],'Color',[0.3686 1.0000 0.4941])% Troll.

hold off
xlim([0 30])
box on
ax.YAxisLocation = 'right';
set(gca,'FontSize',20)
xticklabels('')

nexttile(3)
title('Prairie Potholes and Karsts')
hold on
plot(NaN,NaN,'o--','Color',[0.2118 0.3529 0.4784],'LineWidth',2,'MarkerFaceColor',[0.2118 0.3529 0.4784],'MarkerEdgeColor',[0.2118 0.3529 0.4784])% LV
plot(NaN,NaN,'o-','LineWidth',2,'Color',[0.0588 0.3216 0.1294],'MarkerFaceColor',[0.0588 0.3216 0.1294],'MarkerEdgeColor',[0.0588 0.3216 0.1294])% SE-Deg
plot(NaN,NaN,'o--','LineWidth',2,'MarkerFaceColor',[0.2039 0.6196 0.3137],'MarkerEdgeColor',[0.2039 0.6196 0.3137],'Color',[0.2039 0.6196 0.3137])% CA-DBB
plot(NaN,NaN,'o--','LineWidth',2,'Color',[0.3647 0.7608 0.5765],'MarkerFaceColor',[0.3647 0.7608 0.5765],'MarkerEdgeColor',[0.3647 0.7608 0.5765])% US-ICs
plot(NaN,NaN,'o--','LineWidth',2,'Color',[0.6706 0.7882 0.7176],'MarkerFaceColor',[0.6706 0.7882 0.7176],'MarkerEdgeColor',[0.6706 0.7882 0.7176])% US-Myb
plot(NaN,NaN,'o--','LineWidth',2,'MarkerFaceColor',[0.3686 1.0000 0.4941],'MarkerEdgeColor',[0.3686 1.0000 0.4941],'Color',[0.3686 1.0000 0.4941])% Troll.

plot(binEdges(BG9),meanTA9,'o-','Color',[0.0588 0.4902 0.5216],'MarkerFaceColor',[0.0588 0.4902 0.5216],'MarkerEdgeColor',[0.0588 0.4902 0.5216],'LineWidth',1.5,'MarkerSize',6)% Hogg

for i = 1:length(fn)% Cotton.
    Tableout = newpCO2Tables.(fn{i});
        groupTA11 = discretize(Tableout.AIR_TEMP_C,binEdges);
        [meanTA11,BGZ] = groupsummary(Tableout.CO2_CONC,groupTA11,'mean');
        if length(BGZ(~isnan(BGZ)))>=3
            plot(binEdges(BGZ(~isnan(BGZ))),meanTA11(~isnan(BGZ)),'s-','Color',[0.3059    0.5686    0.6784],'LineWidth',2,'MarkerEdgeColor',[0.3059    0.5686    0.6784],'MarkerFaceColor',[0.3059    0.5686    0.6784])
        end
end

plot(NaN,NaN,'o-k','MarkerFaceColor',[0.7373 0.8196 0.2118],'MarkerEdgeColor',[0.7373 0.8196 0.2118],'LineWidth',1.5,'MarkerSize',6)% SE-Deg
plot(NaN,NaN,'o-','Color',[0.6627 0.9882 0.0118],'MarkerFaceColor',[0.6627 0.9882 0.0118],'MarkerEdgeColor','k','LineWidth',1.5,'MarkerSize',6)% US-Los
plot(NaN,NaN,'o-','Color',[0.3333 0.4902 0.0235],'MarkerFaceColor',[0.3333 0.4902 0.0235],'MarkerEdgeColor',[0.3333 0.4902 0.0235],'LineWidth',1.5,'MarkerSize',6)% US-ALQ (PW)
plot(NaN,NaN,'o--','LineWidth',2,'MarkerFaceColor',[0.3333 0.4902 0.0235],'MarkerEdgeColor','k','Color',[0.3333 0.4902 0.0235])% US-ALQ SW

plot(NaN,NaN,'o-','Color',[0.1451 0.2118 0.0078],'MarkerFaceColor',[0.1451 0.2118 0.0078],'MarkerEdgeColor',[0.1451 0.2118 0.0078],'LineWidth',1.5,'MarkerSize',6)% US-OWC

plot(NaN,NaN,'o--','LineWidth',2,'Color',[0.5020 0.3608 0.5882],'MarkerFaceColor',[0.5020 0.3608 0.5882],'MarkerEdgeColor',[0.5020 0.3608 0.5882])% PHM
plot(NaN,NaN,'o-','LineWidth',2,'Color',[0.7882 0.6000 1.0000],'MarkerFaceColor',[0.7882 0.6000 1.0000],'MarkerEdgeColor',[0.7882 0.6000 1.0000])% US-EDN

hold off
xlim([0 30])
box on
set(gca,'FontSize',20)
ylabel('CO_2({\itaq}) (ppmv)')
xticklabels('')

lg  = legend({'LV','SE-Deg','CA-DBB','US-ICs','US-Myb','Troll.','MBPPW1','Cotton.','','SE-Deg (PW)','US-Los (PW)','US-ALQ (PW)','US-ALQ','US-OWC (PW)','PIE','US-EDN'},...
    'FontSize',20,'Orientation','Horizontal','NumColumns',4);
lg.Layout.Tile = 6;

ax2 = nexttile(4);
title('Porewater')
hold on
plot(binEdges(BGC(~isnan(BGC))),meanTAC(~isnan(BGC)),'o-k','MarkerFaceColor',[0.7373 0.8196 0.2118],'MarkerEdgeColor',[0.7373 0.8196 0.2118],'LineWidth',1.5,'MarkerSize',6)% SE-Deg
plot(binEdges(BG4(~isnan(BG4))),meanTA4(~isnan(BG4)),'o-','Color',[0.6627 0.9882 0.0118],'MarkerFaceColor',[0.6627 0.9882 0.0118],'MarkerEdgeColor','k','LineWidth',1.5,'MarkerSize',6)% US-Los
plot(binEdges(BG(~isnan(meanTA)&~isnan(BG))),meanTA(~isnan(meanTA)&~isnan(BG)),'o-','Color',[0.3333 0.4902 0.0235],'MarkerFaceColor',[0.3333 0.4902 0.0235],'MarkerEdgeColor',[0.3333 0.4902 0.0235],'LineWidth',1.5,'MarkerSize',6)% US-ALQ
plot(binEdges(BG8(~isnan(BG8))),meanTA8(~isnan(BG8)),'o-','Color',[0.1451 0.2118 0.0078],'MarkerFaceColor',[0.1451 0.2118 0.0078],'MarkerEdgeColor',[0.1451 0.2118 0.0078],'LineWidth',1.5,'MarkerSize',6)% US-OWC
hold off
xlim([0 30])
box on
ax2.YAxisLocation = 'right';
set(gca,'FontSize',20)
xlabel('T_{air} (\circC)')

nexttile(5)
hold on
plot(binEdges(BGM(~isnan(BGM))),meanTAM(~isnan(BGM)),'o--','LineWidth',2,'Color',[0.5020 0.3608 0.5882],'MarkerFaceColor',[0.5020 0.3608 0.5882],'MarkerEdgeColor',[0.5020 0.3608 0.5882])% PHM
plot(binEdges(BGX(~isnan(BGX))),meanTAX(~isnan(BGX)),'o-','LineWidth',2,'Color',[0.7882 0.6000 1.0000],'MarkerFaceColor',[0.7882 0.6000 1.0000],'MarkerEdgeColor',[0.7882 0.6000 1.0000])% US-EDN
hold off
xlabel('T_{air} (\circC)')
title('Tidal')
xlim([0 30])
box on
set(gca,'FontSize',20)
ylim([0 5000])

t.TileSpacing = 'compact';
t.Padding = 'compact';

%% WTD and pCO2
binEdges = -0.4:0.2:1.2;

t = tiledlayout(2,2);
title(t,'P_{CO2} vs. WTD','FontSize',20)
ylabel(t,'CO_2({\itaq}) (ppmv)','FontSize',20)
xlabel(t,'WTD (m)','FontSize',20)

nexttile(1);
title('Fens, Bogs, and Marshes')
hold on
SEDeg = SEDeg((ismember(month(SEDeg.DateTime),4:10)),:);
groupTA3 = discretize(SEDeg.WTD,binEdges);
[meanTA3,BG3] = groupsummary(SEDeg.pCO2,groupTA3,'mean');
[R,P] = corrcoef(binEdges(BG3(~isnan(BG3)&~isnan(meanTA3))),meanTA3(~isnan(BG3)&~isnan(meanTA3)));
disp([R,P])
plot(binEdges(BG3(~isnan(BG3)&~isnan(meanTA3))),meanTA3(~isnan(BG3)&~isnan(meanTA3)),'o-','LineWidth',2,'Color',[0.0588 0.3216 0.1294],'MarkerFaceColor',[0.0588 0.3216 0.1294],'MarkerEdgeColor',[0.0588 0.3216 0.1294])% SE-Deg

groupTA2 = discretize(CADBB.WTH*0.001,binEdges);
[meanTA2,BG2] = groupsummary(CADBB.pCO2,groupTA2,'mean');
[R,P] = corrcoef(binEdges(BG2(~isnan(BG2))),meanTA2(~isnan(BG2)));
disp([R,P])
plot(binEdges(BG2(~isnan(meanTA2)&~isnan(BG2))),meanTA2(~isnan(meanTA2)&~isnan(BG2)),'o--','LineWidth',2,'MarkerFaceColor',[0.2039 0.6196 0.3137],'MarkerEdgeColor',[0.2039 0.6196 0.3137],'Color',[0.2039 0.6196 0.3137])% CA-DBB

groupTA8 = discretize(SouthPennines.WTDm,binEdges);
[meanTA8,BG8] = groupsummary(SouthPennines.CO2ppmHeadpace,groupTA8,'mean');
[R,P] = corrcoef(binEdges(BG8(~isnan(BG8))),meanTA8(~isnan(BG8)));
disp([R,P])
plot(binEdges(BG8(~isnan(BG8)&~isnan(meanTA8))),meanTA8(~isnan(BG8)&~isnan(meanTA8)),'o--','LineWidth',2,'Color',[0.4980 0.7608 0.4314],'MarkerFaceColor',[0.4980 0.7608 0.4314],'MarkerEdgeColor',[0.4980 0.7608 0.4314])% South Penn.

groupTAB = discretize(USMyb2.WTD,binEdges);
[meanTAB,BGB] = groupsummary(USMyb2.pCO2,groupTAB,'mean');
[R,P] = corrcoef(binEdges(BGB(~isnan(BGB))),meanTAB(~isnan(BGB)));
disp([R,P])
plot(binEdges(BGB(~isnan(BGB)&~isnan(meanTAB))),meanTAB(~isnan(BGB)&~isnan(meanTAB)),'o--','LineWidth',2,'Color',[0.6706 0.7882 0.7176],'MarkerFaceColor',[0.6706 0.7882 0.7176],'MarkerEdgeColor',[0.6706 0.7882 0.7176])% US-Myb

hold off
xlim([-0.4 1])
box on
set(gca,'FontSize',20)
xticklabels('')

ax = nexttile(2);
title('Prairie Potholes and Karsts')
hold on
plot(NaN,NaN,'o-k','MarkerFaceColor',[0.7373 0.8196 0.2118],'MarkerEdgeColor',[0.7373 0.8196 0.2118],'LineWidth',1.5,'MarkerSize',6)% SE-Deg
plot(NaN,NaN,'o-','LineWidth',2,'Color',[0.0588 0.3216 0.1294],'MarkerFaceColor',[0.0588 0.3216 0.1294],'MarkerEdgeColor',[0.0588 0.3216 0.1294])% SE-Deg
plot(NaN,NaN,'o--','LineWidth',2,'MarkerFaceColor',[0.2039 0.6196 0.3137],'MarkerEdgeColor',[0.2039 0.6196 0.3137],'Color',[0.2039 0.6196 0.3137])% CA-DBB
plot(NaN,NaN,'o--','LineWidth',2,'Color',[0.4980 0.7608 0.4314],'MarkerFaceColor',[0.4980 0.7608 0.4314],'MarkerEdgeColor',[0.4980 0.7608 0.4314])% South Penn.
plot(NaN,NaN,'o--','LineWidth',2,'Color',[0.6706 0.7882 0.7176],'MarkerFaceColor',[0.6706 0.7882 0.7176],'MarkerEdgeColor',[0.6706 0.7882 0.7176])% US-Myb

groupTA1 = discretize(WardBigCypressCO2fluxdataS1.CorrectedStagem,binEdges);
[meanTA1,BG1] = groupsummary(WardBigCypressCO2fluxdataS1.pCO2ppm,groupTA1,'mean');% Big Cypress
[R,P] = corrcoef(binEdges(BG1(~isnan(BG1))),meanTA1(~isnan(BG1)));
disp([R,P])
plot(binEdges(BG1(~isnan(BG1))),meanTA1(~isnan(BG1)),'o--c','LineWidth',2,'MarkerFaceColor','c','MarkerEdgeColor','c')

groupTA = discretize(0.001*MBPPW2.WTD,binEdges);
[meanTA,BG] = groupsummary(MBPPW2.pCO2,groupTA,'mean'); %MBPPW2
[R,P] = corrcoef(binEdges(BG(~isnan(BG))),meanTA(~isnan(BG)));
disp([R,P])
plot(binEdges(BG(~isnan(BG))),meanTA(~isnan(BG)),'o-k','MarkerFaceColor',[0.5961 0.8706 0.8902],'MarkerEdgeColor','k','LineWidth',1.5,'MarkerSize',6)% Young

plot(NaN,NaN,'o-k','MarkerFaceColor',[0.7882 0.6000 1.0000],'MarkerEdgeColor',[0.7882 0.6000 1.0000],'LineWidth',1.5,'MarkerSize',6)% US-EDN

for i = 1:length(fn)
    Tableout = newpCO2Tables.(fn{i});% Cotton.
    groupTA11 = discretize((Tableout.WATER_DEPTH_CM./100),binEdges);
    [meanTA11,BGZ] = groupsummary(Tableout.CO2_CONC,groupTA11,'mean');
    if length(BGZ(~isnan(BGZ)))>=0.4
       [R,P] = corrcoef(binEdges(BGZ(~isnan(BGZ)&~isnan(meanTA11))),meanTA11(~isnan(BGZ)&~isnan(meanTA11)));
       disp([R,P])
       plot(binEdges(BGZ(~isnan(BGZ))),meanTA11(~isnan(BGZ)),'s-','Color',[0.3059    0.5686    0.6784],'LineWidth',2,'MarkerEdgeColor',[0.3059    0.5686    0.6784],'MarkerFaceColor',[0.3059    0.5686    0.6784])
    end
end

hold off
xticklabels('')

xlim([-0.4 1])
box on
set(gca,'FontSize',20)
ax.YAxisLocation = 'right';

lg  = legend({'SE-Deg (PW)','SE-Deg','CA-DBB','South Penn.','US-Myb','Big Cypress','MBPPW2','US-EDN','Cotton.','','','','','','',''},...
    'FontSize',20,'Orientation','Horizontal');
lg.Layout.Tile = 'south';

nexttile(3);
title('Porewater')
hold on
Deg3 = Deg3(ismember(month(Deg3.DateTime),4:10),:);% do not include winter bc water level sensor is trapped in ice
groupTAC = discretize(Deg3.WTD,binEdges);
[meanTAC,BGC] = groupsummary(Deg3.PorewaterCO2,groupTAC,'mean');
[R,P] = corrcoef(binEdges(BGC(~isnan(BGC))),meanTAC(~isnan(BGC)));
disp([R,P])
plot(binEdges(BGC(~isnan(BGC))),meanTAC(~isnan(BGC)),'o-k','MarkerFaceColor',[0.7373 0.8196 0.2118],'MarkerEdgeColor',[0.7373 0.8196 0.2118],'LineWidth',1.5,'MarkerSize',6)% SE-Deg
hold off
xlim([-0.4 1])
box on
set(gca,'FontSize',20)
xlabel('WTD (m)')

ax2 = nexttile(4);
groupTAX = discretize(USEDN.WTD,binEdges);
[meanTAX,BGX] = groupsummary(USEDN.CO2,groupTAX,'mean');
[R,P] = corrcoef(binEdges(BGX(~isnan(BGX))),meanTAX(~isnan(BGX)));
disp([R,P])
plot(binEdges(BGX(~isnan(BGX))),meanTAX(~isnan(BGX)),'o-k','MarkerFaceColor',[0.7882 0.6000 1.0000],'MarkerEdgeColor',[0.7882 0.6000 1.0000],'LineWidth',1.5,'MarkerSize',6)% US-EDN
title('Tidal')
box on
xlim([-0.4 1])
set(gca,'FontSize',20)
xlabel('WTD (m)')

ax2.YAxisLocation = 'right';
t.TileSpacing = 'compact';
t.Padding = 'compact';
