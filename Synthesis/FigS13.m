clear;clc;close all;
load USEDN.mat
%% Convert PCO2 to PCO2BIO
b = 0; % insert calibrated sensitivity parameter
USEDN.CO2 = USEDN.CO2.*(exp(b*(mean(USEDN.Wtemp,'omitnan')-USEDN.Wtemp)));

dailypCO2_mn = [];
dailypCO2_mx = [];
dailydis = [];
dailyGP = [];
dailyRECO = [];

    dailypCO2_min = accumarray(day(USEDN.DateTime,'dayofyear'),USEDN.CO2,[],@min,NaN);
    dailypCO2_max = accumarray(day(USEDN.DateTime,'dayofyear'),USEDN.CO2,[],@max,NaN);
    dailydischarge = accumarray(day(USEDN.DateTime,'dayofyear'),USEDN.water_discharge_m3_2,[],@(x)mean(x,'omitnan'),0);
    dailypCO2_mn = [dailypCO2_mn; dailypCO2_min];
    dailypCO2_mx = [dailypCO2_mx; dailypCO2_max];
    dailydis = [dailydis; dailydischarge];
    dailyGPP = accumarray(day(USEDN.DateTime,'dayofyear'),USEDN.GPP_DT_uStar,[],@(x)mean(x,'omitnan'),0);
    dailyGP = [dailyGP; dailyGPP];
    dailyReco = accumarray(day(USEDN.DateTime,'dayofyear'),USEDN.Reco_DT_uStar,[],@(x)mean(x,'omitnan'),0);
    dailyRECO = [dailyRECO; dailyReco];

dailypCO2_bio = dailypCO2_mx-dailypCO2_mn;
DateTime = USEDN.DateTime(1:48:end);
DateTime = dateshift(DateTime, 'start', 'day');
dailydis = dailydis(dailydis~=0);
%% comparison of hourly-scale data
% Calculate correlation coefficient between PCO2 and GPP for each day
% reshape PCO2 and GPP into daily rows, calculate corrcoef
r = rem(length(USEDN.CO2),48);

dailypCO2 = reshape(USEDN.CO2(1:(end-r)),48,[]);
dailyGPP = reshape(USEDN.GPP_DT_uStar(1:(end-r)),48,[]);
R_all = zeros(1,105);
P_all = zeros(1,105);
count = 1;

for i = 1:width(dailyGPP)
    A = dailypCO2(:,i);
    B = dailyGPP(:,i);
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

%% Regimes

d = dailydis(1:end-2);
P = prctile(d,[25 50 75]);
flows = {'no','low','medium','high'};
binnedData = discretize(d,[min(dailydis) P(1) P(2) P(3) max(dailydis)],'categorical',flows);
hr = hour(DateTime(1:end-1));

highflowR = R_all(binnedData=='high');
highflowP = P_all(binnedData=='high');
h = highflowR(~isnan(highflowR)&(highflowP<0.05));
mean(h)
width(h)
mean(highflowP(~isnan(highflowR)&(highflowP<0.05)))

medflowR = R_all(binnedData=='medium');
medflowP = P_all(binnedData=='medium');
m = medflowR(~isnan(medflowR)&(medflowP<0.05));
mean(m)
width(m)
mean(medflowP(~isnan(medflowR)&(medflowP<0.05)))

lowflowR = R_all(binnedData=='low');
lowflowP = P_all(binnedData=='low');
l = lowflowR(~isnan(lowflowR)&(lowflowP<0.05));
mean(l)
width(l)
mean(lowflowP(~isnan(lowflowR)&(lowflowP<0.05)))

noflowR = R_all(binnedData=='no');
noflowP = P_all(binnedData=='no');
n = noflowR(~isnan(noflowR)&(noflowP<0.05));
mean(n)
width(n)
mean(noflowP(~isnan(noflowR)&(noflowP<0.05)))
 
x = d(P_all<0.05);
y = R_all(P_all<0.05);
x(isnan(y))= nan;
y(isnan(x))=nan;
[R,P]=corrcoef(x(~isnan(x)),(y(~isnan(y))));
%% Plotting

figure()
plot(d(P_all<0.05),R_all(P_all<0.05),'.k','MarkerSize',20)
xlabel('discharge (m^3s^{-1})')
ylabel('corr. coeff.')
title('US-EDN Discharge vs. GPP & P_{CO2} Correlation')
set(gca,'FontSize',15)

