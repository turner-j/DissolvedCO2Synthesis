% Import daily average data
clear; clc;close all;
addpath(genpath(pwd))
load USEDNPCO2_daily.mat
load USALQStreamDaily.mat
USALQ = readtable("E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls",'FileType','spreadsheet','Sheet','US-ALQ_DT');
USLos = readtable("E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls",'FileType','spreadsheet','Sheet','US-Los_DT');
CADBB2 = readtable("E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls",'FileType','spreadsheet','Sheet','CA-DBB2');
SEDeg = readtable("E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls",'FileType','spreadsheet','Sheet','SEDeg_DT');% not the interpolated version
SEDeg2 = readtable("E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls",'FileType','spreadsheet','Sheet','SEDeg_DTinterp');
SEDegPW = readtable("E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls",'FileType','spreadsheet','Sheet','PorewaterSEDeg_DT');
USGCE = readtable("E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls",'FileType','spreadsheet','Sheet','US-GCE_DT');
LochVale = readtable("E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls",'FileType','spreadsheet','Sheet','LochVale');
load('LochValepCO2.mat')
USOWC = readtable("E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls",'FileType','spreadsheet','Sheet','US-OWC_DT');
USEvM = readtable("E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls",'FileType','spreadsheet','Sheet','USEvM_DT');
MBPPW1_interp = readtable("E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls",'FileType','spreadsheet','Sheet','MBPPW1FCO2_interp');
MBPPW2_interp = readtable("E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls",'FileType','spreadsheet','Sheet','MBPPW2FCO2_interp');
MBPPW1 = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','MBPPW1FCO2');
MBPPW2 = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','MBPPW2FCO2');

addpath('E:\SanDiskSecureAccess\all MATLAB files\Residence paper\YKD')
YKDUB = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','YKDUB');
YKDB = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','YKDB');
USICs = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','USICs_DT');

PorewaterSEDeg = readtable("E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls",'FileType','spreadsheet','Sheet','PorewaterSEDeg_DT');
load USMyb.mat;USMyb2 = USMyb2(145:end,:);
Troll = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','FileType','spreadsheet','Sheet','Troll_interp');
USICs_interp = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','FileType','spreadsheet','Sheet','USICs_DTinterp');
BZF = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','FileType','spreadsheet','Sheet','APEX_interp');
BZF2 = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','FileType','spreadsheet','Sheet','APEX');

PLM = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','PLM2');
PHM = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','PHM2');
USHB1 = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','USHB1');

LochValeGPP = 0.063675;

USOWC = USOWC(6:end,:);
OWCsoilCO2 = USOWC(:,12:end);
USOWC.pCO2 = mean(table2array(OWCsoilCO2),2,'omitnan');% take the mean of each row
load USEDN.mat
%%
load SouthPennines.mat
run ImportYKDpCO2.m
run ImportBigCypresspCO2.m

YKDpCO2burned_fen = ZolkosetalYKDaquaticsfinal((ZolkosetalYKDaquaticsfinal.LandscapeCategory==('fen'))&(ZolkosetalYKDaquaticsfinal.WaterType=='surface')&(ZolkosetalYKDaquaticsfinal.Burn=='Burned'),:);
YKDpCO2burned_bog = ZolkosetalYKDaquaticsfinal((ZolkosetalYKDaquaticsfinal.LandscapeCategory==('pond')|(ZolkosetalYKDaquaticsfinal.LandscapeCategory=='peat plateau'))&(ZolkosetalYKDaquaticsfinal.WaterType=='surface')&(ZolkosetalYKDaquaticsfinal.Burn=='Burned'),:);

YKDpCO2burned_fenPW = ZolkosetalYKDaquaticsfinal((ZolkosetalYKDaquaticsfinal.LandscapeCategory==('fen'))&(ZolkosetalYKDaquaticsfinal.WaterType=='pore')&(ZolkosetalYKDaquaticsfinal.Burn=='Burned'),:);
YKDpCO2burned_bogPW = ZolkosetalYKDaquaticsfinal((ZolkosetalYKDaquaticsfinal.LandscapeCategory==('pond')|(ZolkosetalYKDaquaticsfinal.LandscapeCategory=='peat plateau'))&(ZolkosetalYKDaquaticsfinal.WaterType=='pore')&(ZolkosetalYKDaquaticsfinal.Burn=='Burned'),:);

YKDpCO2unburned_fen = ZolkosetalYKDaquaticsfinal((ZolkosetalYKDaquaticsfinal.LandscapeCategory==('fen'))&(ZolkosetalYKDaquaticsfinal.WaterType=='surface')&(ZolkosetalYKDaquaticsfinal.Burn=='Unburned'),:);
YKDpCO2unburned_bog = ZolkosetalYKDaquaticsfinal((ZolkosetalYKDaquaticsfinal.LandscapeCategory==('pond')|(ZolkosetalYKDaquaticsfinal.LandscapeCategory=='peat plateau'))&(ZolkosetalYKDaquaticsfinal.WaterType=='surface')&(ZolkosetalYKDaquaticsfinal.Burn=='Unburned'),:);

YKDpCO2unburned_fenPW = ZolkosetalYKDaquaticsfinal((ZolkosetalYKDaquaticsfinal.LandscapeCategory==('fen'))&(ZolkosetalYKDaquaticsfinal.WaterType=='pore')&(ZolkosetalYKDaquaticsfinal.Burn=='Unburned'),:);
YKDpCO2unburned_bogPW = ZolkosetalYKDaquaticsfinal((ZolkosetalYKDaquaticsfinal.LandscapeCategory==('pond')|(ZolkosetalYKDaquaticsfinal.LandscapeCategory=='peat plateau'))&(ZolkosetalYKDaquaticsfinal.WaterType=='pore')&(ZolkosetalYKDaquaticsfinal.Burn=='Unburned'),:);

USALQ.pCO2(USALQ.pCO2<=416) = NaN;% 2021

WetlandType = categorical();
for i = 1:size(USALQStreamDaily,1)
    WetlandType(i) =  categorical({'Fen'});% hourly measurements
end
USALQStreamDaily.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(CADBB2)
    WetlandType(i) =  categorical({'Bog'});
end
CADBB2.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(SEDeg,1)
    WetlandType(i) =  categorical({'Fen'});
end
SEDeg.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(USGCE,1)
    WetlandType(i) =  categorical({'Tidal'});
end
USGCE.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(LochVale,1)
    WetlandType(i) =  categorical({'Alpine'});
end
LochVale.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(USEvM,1)
    WetlandType(i) =  categorical({'Tidal'});
end
USEvM.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(MBPPW1.pCO2,1)
    WetlandType(i) =  categorical({'Prairie Pothole or Karst'});
end
MBPPW1.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(MBPPW2,1)
    WetlandType(i) =  categorical({'Prairie Pothole or Karst'});
end
MBPPW2.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(USICs,1)
    WetlandType(i) =  categorical({'Fen'});
end
USICs.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(USMyb2,1)
    WetlandType(i) =  categorical({'Marsh'});
end
USMyb2.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(Troll,1)
    WetlandType(i) =  categorical({'Fen'});
end
Troll.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(PLM,1)
    WetlandType(i) =  categorical({'Tidal'});
end
PLM.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(USHB1,1)
    WetlandType(i) =  categorical({'Tidal'});
end
USHB1.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(USEDNPCO2DAY,1)
    WetlandType(i) =  categorical({'Tidal'});
end
USEDNPCO2DAY.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(YKDpCO2burned_fen,1)
    WetlandType(i) =  categorical({'Fen'});
end
YKDpCO2burned_fen.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(YKDpCO2unburned_fen,1)
    WetlandType(i) =  categorical({'Fen'});
end
YKDpCO2unburned_fen.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(YKDpCO2burned_bog,1)
    WetlandType(i) =  categorical({'Bog'});
end
YKDpCO2burned_bog.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(YKDpCO2unburned_bog,1)
    WetlandType(i) =  categorical({'Bog'});
end
YKDpCO2unburned_bog.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(SouthPennines,1)
    WetlandType(i) =  categorical({'Bog'});
end
SouthPennines.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(WardBigCypressCO2fluxdataS1,1)
    WetlandType(i) =  categorical({'Prairie Pothole or Karst'});
end
WardBigCypressCO2fluxdataS1.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(YKDUB,1)
    WetlandType(i) =  categorical({'Fen'});
end
YKDUB.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(YKDB,1)
    WetlandType(i) =  categorical({'Fen'});
end
YKDB.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(USEDNPCO2DAY,1)
    WetlandType(i) =  categorical({'Tidal'});
end
USEDNPCO2DAY.WetlandType = WetlandType';

y2 = [USALQStreamDaily.pCO2_ppm;...
    CADBB2.pCO2;SEDeg.pCO2;...
    USGCE.pCO2;LochVale.pCO2;...
    USEvM.pCO2;MBPPW1.pCO2;MBPPW2.pCO2;USICs.pCO2;...
    USMyb2.pCO2;Troll.pCO2;...
    PLM.pCO2;USHB1.pCO2;USEDNPCO2DAY.CO2;...
    YKDpCO2burned_fen.pCO2;YKDpCO2unburned_fen.pCO2;...
    YKDpCO2burned_bog.pCO2;YKDpCO2unburned_bog.pCO2;...
    SouthPennines.CO2ppmHeadpace;WardBigCypressCO2fluxdataS1.pCO2ppm;...
    YKDUB.pCO2;YKDB.pCO2;USEDNPCO2DAY.CO2];

g2 = [USALQStreamDaily.WetlandType;
    CADBB2.WetlandType;SEDeg.WetlandType;
    USGCE.WetlandType;LochVale.WetlandType;
    USEvM.WetlandType;MBPPW1.WetlandType;MBPPW2.WetlandType;USICs.WetlandType;
    USMyb2.WetlandType;Troll.WetlandType;
    PLM.WetlandType;USHB1.WetlandType;USEDNPCO2DAY.WetlandType;
    YKDpCO2burned_fen.WetlandType;YKDpCO2unburned_fen.WetlandType;
    YKDpCO2burned_bog.WetlandType;YKDpCO2unburned_bog.WetlandType;
    SouthPennines.WetlandType;WardBigCypressCO2fluxdataS1.WetlandType;
    YKDUB.WetlandType;YKDB.WetlandType;USEDNPCO2DAY.WetlandType];

y = [BZF2.pCO2;...
    USALQ.pCO2;SEDegPW.PorewaterCO2;USLos.pCO2;...
    YKDpCO2burned_fenPW.pCO2;YKDpCO2burned_bogPW.pCO2;YKDpCO2unburned_fenPW.pCO2;...
    YKDpCO2unburned_bogPW.pCO2;...
    USOWC.pCO2];

% Porewater 
WetlandType = categorical();
for i = 1:size(BZF2,1)
    WetlandType(i) =  categorical({'Fen'});
end
BZF2.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(USALQ,1)
    WetlandType(i) =  categorical({'Fen'});
end
USALQ.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(SEDegPW,1)
    WetlandType(i) =  categorical({'Fen'});
end
SEDegPW.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(USLos,1)
    WetlandType(i) =  categorical({'Fen'});
end
USLos.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(YKDpCO2burned_fenPW,1)
    WetlandType(i) =  categorical({'Fen'});
end
YKDpCO2burned_fenPW.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(YKDpCO2burned_bogPW,1)
    WetlandType(i) =  categorical({'Bog'});
end
YKDpCO2burned_bogPW.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(YKDpCO2unburned_fenPW,1)
    WetlandType(i) =  categorical({'Fen'});
end
YKDpCO2unburned_fenPW.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(YKDpCO2unburned_bogPW,1)
    WetlandType(i) =  categorical({'Bog'});
end
YKDpCO2unburned_bogPW.WetlandType = WetlandType';

WetlandType = categorical();
for i = 1:size(USOWC,1)
    WetlandType(i) =  categorical({'Marsh'});
end
USOWC.WetlandType = WetlandType';

t = tiledlayout(1,3);
nexttile(1,[1 2]);
[means,grps] = grpstats(y2,g2, ...
    ["mean","gname"]);
x = reordercats(categorical(grps),{'Alpine' 'Marsh' 'Tidal' 'Fen' 'Bog' 'Prairie Pothole or Karst'});
b = bar(x,means,'EdgeColor','none');

hold on
[sem,~] = grpstats(y2,g2, ...
    ["sem","gname"]);
er = errorbar(x,means,sem);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 2;

hold off

ylabel('CO_2({\itaq}) (ppmv)')
b.FaceColor = 'flat';
b.CData(1,:) = [0.2118 0.3529 0.4784];% Alpine
b.CData(2,:) = [.8471 0.6706 1.0000];% Marsh
b.CData(3,:) = [0.5020 0.3608 0.5882];% Tidal
b.CData(4,:) = [0.2039 0.6196 0.3137];% Fen
b.CData(5,:) = [0.5451 0.8000 0.3765];% Bog
b.CData(6,:) = [0.5961 0.8706 0.8902];% Prairie Pothole or Karst
set(gca,'TickLength',[0 .01],'FontSize',20)

y = y(~isnan(y));
g = categorical();
for i = 1:length(y)
    g(i) =  categorical({'Porewater'});
end

ax = nexttile(3);
[means,grps] = grpstats(y,g, ...
    ["mean","gname"]);
x = reordercats(categorical(grps),{'Porewater'});
b = bar(x,means,'EdgeColor','none');
ylabel('')
b.FaceColor = 'flat';
b.CData = [0.3333 0.4902 0.0235];% Porewater
set(gca,'FontSize',12)

hold on
[sem,~] = grpstats(y,g, ...
    ["sem","gname"]);
er = errorbar(x,means,sem);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 2;

hold off

ax.YAxisLocation = 'right';
set(gca,'TickLength',[0 .01],'FontSize',20)
xtickangle(45)

t.TileSpacing = 'compact';
t.Padding = 'compact';

[p,tbl,stats]  = anova1(y2,g2);
disp(p)
results = multcompare(stats);