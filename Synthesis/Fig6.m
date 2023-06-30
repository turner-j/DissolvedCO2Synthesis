clear;clc;

addpath 'E:\SanDiskSecureAccess\all MATLAB files\allequash'
addpath 'E:\SanDiskSecureAccess\all MATLAB files\Residence paper\USALQ'
addpath 'C:\Users\lturn\Desktop\all MATLAB files\Residence paper\NewAnalysis_2023'
run ImportALQREddyProcOut2.m
run ALQstreampCO2_2023.m
load('USALQ2021.mat')
%% Create matching datetime variables 
USALQ.NEE_F = ALQREddyProcOut.NEE_uStar_f;
USALQ.GPP_F = ALQREddyProcOut.Reco_DT_uStar;

% preallocate variable outcome for speed
X = zeros(365,38);
count = 1;

for j= 3:width(USALQ)
    USALQ.(j)(USALQ.(j)==-9999) = NaN;
    dailyave = dailyaverage2(USALQ.(j));
    X(:,count) = dailyave(:,1);
    count = count + 1;
end

% add variable names, timestamp
X = array2table(X);
X.Properties.VariableNames = USALQ.Properties.VariableNames;
X.TIMESTAMP_END = USALQ.TIMESTAMP_END(dailyave(:,2).*48);

%% Create daily averages by indexing
yr = year(Allequashstream2023.Date_Time);
Allequashstream2023 = Allequashstream2023(yr==2021,:);

%% Convert PCO2 to PCO2BIO
Allequashstream2023.Discharge_CMS = Allequashstream2023.Discharge_CFS.*0.028316847; % convert cfs to m3/s
b = 0.0290;% fill in with the calibrated sensitivity parameter

% Water-temperature-normalize PCO2
Allequashstream2023.pCO2_ppm_norm = Allequashstream2023.pCO2_ppm.*(exp(b*(mean(Allequashstream2023.Water_Temp_C,'omitnan')-Allequashstream2023.Water_Temp_C)));
dailypCO2_min = accumarray(day(Allequashstream2023.Date_Time,'dayofyear'),Allequashstream2023.pCO2_ppm,[],@min,NaN);
dailypCO2_max = accumarray(day(Allequashstream2023.Date_Time,'dayofyear'),Allequashstream2023.pCO2_ppm,[],@max,NaN);
dailypCO2_bio = dailypCO2_max-dailypCO2_min;

dailydischarge = accumarray(day(Allequashstream2023.Date_Time,'dayofyear'),Allequashstream2023.Discharge_CMS,[],@(x)mean(x,'omitnan'),0);

[~,ia] = unique(day(Allequashstream2023.Date_Time,'dayofyear'));
dailypco2datesALQ = dateshift(Allequashstream2023.Date_Time(ia), 'start', 'day');

X.SWpCO2 = nan(size(X,1),1);X.discharge = nan(size(X,1),1);
ix = (ismember(dailypco2datesALQ,X.TIMESTAMP_END));
X.SWpCO2(ix) = dailypCO2_bio(ix(ix~=0));
X.discharge(ix) = dailydischarge(ix(ix~=0));

USALQ.SWpCO2 = nan(size(USALQ,1),1);USALQ.discharge = nan(size(USALQ,1),1);
ik = (ismember(USALQ.TIMESTAMP_END,Allequashstream2023.Date_Time));
ik2 = (ismember(Allequashstream2023.Date_Time,USALQ.TIMESTAMP_END));
USALQ.SWpCO2(ik) = Allequashstream2023.pCO2_ppm_norm(ik2);
USALQ.discharge(ik) = Allequashstream2023.Discharge_CMS(ik2);

%% Plotting comparison of hourly-scale data
% Calculate correlation coefficient between PCO2 and GPP for each day
% reshape PCO2 and GPP into daily rows, calculate corrcoef

dailypCO2 = reshape(USALQ.SWpCO2,48,[]);
dailyGPP = reshape(USALQ.GPP_F,48,[]);
R_all = zeros(1,365);
P_all = zeros(1,365);
count = 1;

for i = 1:width(dailyGPP)
    A = dailypCO2(:,i);
    B = dailyGPP(:,i);
    [R,P] = corrcoef(A(~isnan(A)&~isnan(B)),B(~isnan(B)&~isnan(A)));
    R_all(count) = R(2);
    P_all(count) = P(2);
    count = count + 1;
end

R_all(R_all<0) = NaN;

z = tiledlayout(2,3,'TileSpacing','Compact');

nexttile(1)
plot(X.TIMESTAMP_END(~isnan(X.discharge)),X.discharge(~isnan(X.discharge)),'-k','LineWidth',2)
xlim([X.TIMESTAMP_END(1) X.TIMESTAMP_END(274)])
ylabel('discharge (m^{3}s^{-1})')
title('US-ALQ')
xticklabels({''})
set(gca,'FontSize',22)
ylim([0 0.5])

nexttile(4)
hold on
t = USALQ.TIMESTAMP_END(1:48:end);
scatter(t,R_all,100,'MarkerFaceColor','w','MarkerEdgeColor','k')
scatter(t(P_all<0.05),R_all(P_all<0.05),100,'MarkerFaceColor',[0.3333 0.4902 0.0235],'MarkerEdgeColor','k')
hold off
xlim([X.TIMESTAMP_END(1) X.TIMESTAMP_END(274)])
ylabel('corr. coeff.')
box on
set(gca,'FontSize',22)
xtickangle(45)
%% US-EDN
load USEDN.mat
%% Convert PCO2 to PCO2BIO
b = 0; % insert calibrated sensitivity parameter
USEDN.CO2 = USEDN.CO2.*(exp(b*(mean(USEDN.Wtemp,'omitnan')-USEDN.Wtemp)));

dailydischarge = accumarray(day(USEDN.DateTime,'dayofyear'),USEDN.water_discharge_m3_2,[],@(x)mean(x,'omitnan'),NaN);
DateTime = USEDN.DateTime(1:48:end);
DateTime = dateshift(DateTime, 'start', 'day');
dailydischarge = dailydischarge(~isnan(dailydischarge));
%% Plotting comparison of hourly-scale data
% Calculate correlation coefficient between PCO2 and GPP for each day
% reshape PCO2 and GPP into daily rows, calculate corrcoef
r = rem(length(USEDN.CO2),48);

dailypCO2 = reshape(USEDN.CO2(1:(end-r)),48,[]);
dailyReco = reshape(USEDN.Reco_DT_uStar(1:(end-r)),48,[]);

clear R_all P_all
R_all = zeros(1,105);
P_all = zeros(1,105);
count = 1;

for i = 1:width(dailyReco)
    A = dailypCO2(:,i);
    B = dailyReco(:,i);
    [R,P] = corrcoef(A(~isnan(A)&~isnan(B)),B(~isnan(B)&~isnan(A)));
    if length(R)>1
    R_all(count) = R(2);
    P_all(count) = P(2);
    else 
    R_all(count) = NaN;
    P_all(count) = NaN;
    end
    count = count + 1;
end

R_all(R_all<0) = NaN;

nexttile(2)
plot(DateTime,dailydischarge,'-k','LineWidth',2)
xlim([DateTime(1) DateTime(end)])
title('US-EDN')
xticklabels({''})
set(gca,'FontSize',20)

nexttile(5)
hold on
DTime = DateTime(1:end-1);
scatter(DTime,R_all,100,'MarkerFaceColor','w','MarkerEdgeColor','k')
scatter(DTime(P_all<0.05),R_all(P_all<0.05),100,'MarkerFaceColor',[0.7882 0.6000 1.0000],'MarkerEdgeColor','k')
hold off
xlim([DateTime(1) DateTime(end)])
box on
ylim([0 1])
set(gca,'FontSize',20)
%% SE-Deg
cd('E:\SanDiskSecureAccess\all MATLAB files\Residence paper\SE-Deg')
DEG_Discharge = readtable("E:\SanDiskSecureAccess\all MATLAB files\Residence paper\SE-Deg\StreamData.csv",'ReadVariableNames', false,'HeaderLines',28);
run ImportHourlyWaterTable.m

pCO2_Deg = readtable('C:\Users\lturn\Downloads\20221103-Degero-mire.xlsx');
run ImportSEDegREddyProcOut3.m
Degoriginal = readtable('C:\Users\lturn\Downloads\FLX_SE-Deg_FLUXNET-CH4_2014-2018_1-1\FLX_SE-Deg_FLUXNET-CH4_2014-2018_1-1\FLX_SE-Deg_FLUXNET-CH4_DD_2014-2018_1-1.csv');
Degoriginal.TIMESTAMP = datetime(Degoriginal.TIMESTAMP,'ConvertFrom','yyyyMMdd','Format','yyyyMMdd');

% Load porewater data
run ImportDailyPorewaterCO2_ppm.m
%% Index for days with concurrent data
degeroRdatabase2020 = degeroRdatabase2020(~isnat(degeroRdatabase2020.DATE_TIME),:);
[C,ia] = (unique(degeroRdatabase2020.DATE_TIME));
degeroRdatabase2020 = degeroRdatabase2020(ia,:);

SEDegREddyProcOut2.pCO2 = NaN(size(SEDegREddyProcOut2,1),1);
ix =find(ismember(SEDegREddyProcOut2.DateTime,degeroRdatabase2020.DATE_TIME));

degeroRdatabase2020.maxpCO2 = max(table2array(degeroRdatabase2020(:,2:end)),[],2);
SEDegREddyProcOut2.pCO2(ix)=degeroRdatabase2020.maxpCO2;

SEDegREddyProcOut2.Twater = NaN(size(SEDegREddyProcOut2,1),1);
ixx = find(ismember(SEDegREddyProcOut2.DateTime,HourlyPorewaterCO2Temp.DATE_TIME));
ind = find(ismember(HourlyPorewaterCO2Temp.DATE_TIME,SEDegREddyProcOut2.DateTime(ixx)));
SEDegREddyProcOut2.Twater(ixx) = HourlyPorewaterCO2Temp.Temp_75(ind,:);
SEDegREddyProcOut2 = SEDegREddyProcOut2(7915:70129,:);
%% Convert PCO2 to PCO2BIO
DEG_Discharge.DischargeCMS = DEG_Discharge.Var3.*0.001; % convert L/s to m3/s
SEDegREddyProcOut2.DischargeCMS = NaN(size(SEDegREddyProcOut2,1),1);

ind = find(ismember(DEG_Discharge.Var2,SEDegREddyProcOut2.DateTime));
ikx =find(ismember(SEDegREddyProcOut2.DateTime,DEG_Discharge.Var2));
SEDegREddyProcOut2.DischargeCMS(ikx) = DEG_Discharge.DischargeCMS(ind);

%% Water-temperature-normalize PCO2
b = 0.0230; % insert calibrated sensitivity parameter

SEDegREddyProcOut2.pCO2_ppm_norm = SEDegREddyProcOut2.pCO2.*(exp(b*(mean(SEDegREddyProcOut2.Twater,'omitnan')-SEDegREddyProcOut2.Twater)));
yr = year(SEDegREddyProcOut2.DateTime);
dailypCO2_mn = [];
dailypCO2_mx = [];
dailydis = [];
dailyGP = [];

for i = yr(1):yr(end)
    ind = find(ismember(yr,i));
    Deg = SEDegREddyProcOut2(ind,:);
    dailypCO2_min = accumarray(day(Deg.DateTime,'dayofyear'),Deg.pCO2_ppm_norm,[],@min,NaN);
    dailypCO2_max = accumarray(day(Deg.DateTime,'dayofyear'),Deg.pCO2_ppm_norm,[],@max,NaN);
    dailydischarge = accumarray(day(Deg.DateTime,'dayofyear'),Deg.DischargeCMS,[],@(x)mean(x,'omitnan'),0);
    dailypCO2_mn = [dailypCO2_mn; dailypCO2_min];
    dailypCO2_mx = [dailypCO2_mx; dailypCO2_max];
    dailydis = [dailydis; dailydischarge];
    dailyGPP = accumarray(day(Deg.DateTime,'dayofyear'),Deg.Reco_DT_uStar,[],@(x)mean(x,'omitnan'),0);
    dailyGP = [dailyGP; dailyGPP];
end

dailypCO2_bio = dailypCO2_mx-dailypCO2_mn;
DateTime = SEDegREddyProcOut2.DateTime(1:48:end);
DateTime = dateshift(DateTime, 'start', 'day');
dailydis = dailydis(dailydis~=0);
%% Plotting comparison of hourly-scale data
% Calculate correlation coefficient between PCO2 and GPP for each day
% reshape PCO2 and GPP into daily rows, calculate corrcoef
r = rem(length(SEDegREddyProcOut2.pCO2_ppm_norm),48);

dailypCO2 = reshape(SEDegREddyProcOut2.pCO2_ppm_norm(1:(end-r)),48,[]);
dailyGPP = reshape(SEDegREddyProcOut2.GPP_DT_uStar(1:(end-r)),48,[]);

clear R_all P_all
R_all = zeros(1,1296);
P_all = zeros(1,1296);

for i = 1:width(dailyGPP)
    A = dailypCO2(:,i);
    B = dailyGPP(:,i);
    [R,P] = corrcoef(A(~isnan(A)&~isnan(B)),B(~isnan(B)&~isnan(A)));
    if length(R)>1
    R_all(i) = R(2);
    P_all(i) = P(2);
    else 
    R_all(i) = NaN;
    P_all(i) = NaN;
    end
end

R_all(R_all<0) = NaN;

for j = unique(year(DateTime))'
    DT = DateTime(year(DateTime)==j);
    nexttile(3)
    plot(DateTime,dailydis(1:end-1),'-k','LineWidth',2)
    xlim([DT(1) DT(end)])
    title('SE-Deg')
    xticklabels({''})
    ylim([0 0.5])
    set(gca,'FontSize',20)

    nexttile(6)
    hold on
    DTime = DateTime(1:end-1);
    scatter(DTime,R_all,100,'MarkerFaceColor','w','MarkerEdgeColor','k')
    scatter(DTime(P_all<0.05),R_all(P_all<0.05),100,'MarkerFaceColor',[0.7373 0.8196 0.2118],'MarkerEdgeColor','k')
    hold off
    xlim([DT(1) DT(end)])
    box on
    set(gca,'FontSize',20)
    xtickangle(45)
end

