  clear;clc;
load('LochValepCO2.mat')
run ImportLVMetData.m
addpath 'E:\SanDiskSecureAccess\all MATLAB files\allequash'
load("LVDailyData.mat")
load('USGS401723105400000.mat')
LochValeClimate = WY1992to2019LochValeClimate((WY1992to2019LochValeClimate.station_name=='Andrews Creek weather station')&(WY1992to2019LochValeClimate.measurement_height==4),:);

ind = find(ismember(LochValeClimate.timestamp,AndrewsHourlyData.DateTime));

LochValeClimate = LochValeClimate(ind,:);
LochValeClimate.pCO2 = AndrewsHourlyData.pCO2_uatm;
s = string(LochValeClimate.Snow_depth);LochValeClimate.Snow_depth = double(s);

LVDailyData = LVDailyData(LVDailyData.Site_ID==401723105400000,:);

%% Save half-hourly data then take daily averages
i = 1;LochVale = [];
LochValeClimate = LochValeClimate(year(LochValeClimate.timestamp)==2017,:);
for j= [2:8 10]
    LochVale(:,i) = accumarray(day(LochValeClimate.timestamp,'dayofyear'),LochValeClimate.(j),[],@mean);
    i = i + 1;
end

% add variable names, timestamp
LochVale = array2table(LochVale);
LochVale.Properties.VariableNames = LochValeClimate.Properties.VariableNames([2:8 10]);
newtimes = unique(day(LochValeClimate.timestamp,'dayofyear'));
%% Plot pCO2 and FCO2
% Convert mmol/m2/d to umol/m2/s
FCO2 = (Mastetal1998LochValeCO2flux.SaturatedAverageDailyCO2Fluxmmolm2d1.*1000)./86400;
LochValeClimate.Snow_depth(isnan(LochValeClimate.Snow_depth)) = 0;
c = LochVale.Snow_depth;

f = figure();
right_color = [0 0 0];
left_color = [0 0 0];
set(f,'defaultAxesColorOrder',[left_color; right_color]);

[M,I] = max(LochValeClimate.Snow_depth);
monthlysnow = accumarray(month(LochValeClimate.timestamp),LochValeClimate.Snow_depth,[],@mean,NaN);

yyaxis left
hold on
sz = (c+1)*50;
scatter(newtimes,LochVale.pCO2,sz,c,'filled','MarkerEdgeColor','k')
hold off
ylabel('CO_2({\itaq}) (ppmv)')
datetick('x','mmm')
yyaxis right
plot(day(Mastetal1998LochValeCO2flux.Date,'dayofyear'),FCO2,'-k','MarkerSize',12,'LineWidth',2)
ylabel('CO_2 (\mumol CO_2 m^-^2s^-^1)')
legend('P_{CO2} 2017','CO_2 1996','Location','northeast','NumColumns',2)
c = colorbar;
colormap("parula")
colormap(flipud(colormap))
c.Label.String = 'snow depth (m)';
c.Label.FontSize = 20;
ax = gca; 
ax.FontSize = 20;
xlim([1 200])
