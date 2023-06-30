% Clear workspace
close all;clear; clc;
addpath(genpath(pwd))
load USMyb.mat
USMybREddyProcout = USMybREddyProcout(7195:end,:);
USMybREddyProcout.pCO2(USMybREddyProcout.pCO2<405.0)=NaN;
load SEDegPorewater.mat
SEDegREddyProcOut2(1:7914,:) = [];
SEDegREddyProcOut2(2:2:end,:)=[];
LochVale = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','LochVale');
USLos2 = readtable("E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls",'FileType','spreadsheet','Sheet','US-Los');
USALQ2 = readtable("E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls",'FileType','spreadsheet','Sheet','US-ALQ_DT');
USALQ = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','US-ALQ_DT');
ind = find(ismember(USALQ2.TIMESTAMP_END,USALQ.TIMESTAMP_END));
USALQ2 = USALQ2(ind(1):ind(end),:);
ind = find(ismember(USALQ.TIMESTAMP_END,USALQ2.TIMESTAMP_END));
USALQ = USALQ(ind(1):ind(end),:);

USLos2 = USLos2((ismember(day(USLos2.TIMESTAMP_END,'dayofyear'),229:366)&ismember(year(USLos2.TIMESTAMP_END),2020)),:);
USLos2.pCO2(USLos2.pCO2<413)=NaN;
USALQ2.pCO2(USALQ2.pCO2<415)=NaN;
USLos = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','US-Los_DT');
USLos = USLos((ismember(day(USLos.TIMESTAMP_END,'dayofyear'),229:366)&ismember(year(USLos.TIMESTAMP_END),2020)),:);
USLos = USLos(5:132,:);
%% Plotting
f = tiledlayout(2,2);
% Yearly
nexttile(3)

% US-Myb
t = zeros(32,2);
count = 1;
hold on 

yearlypCO2_orig = accumarray(year(USMyb2.DateTime),USMyb2.pCO2,[],@(x)mean(x,'omitnan'),NaN);

for i = 1:10:1000
    y = USMyb2((1:i:end),:);
    yearlypCO2 = accumarray(year(y.DateTime),y.pCO2,[],@(x)mean(x,'omitnan'),NaN);
    if length(yearlypCO2(~isnan(yearlypCO2)))==length(yearlypCO2_orig(~isnan(yearlypCO2_orig)))
        R = corrcoef(yearlypCO2(~isnan(yearlypCO2)&~isnan(yearlypCO2_orig)),yearlypCO2_orig(~isnan(yearlypCO2)&~isnan(yearlypCO2_orig)));
        t(count,1) = 365/i;
        t(count,2) = R(2);
        count = count + 1;
    end
end

z = t;

plot(t(:,1),t(:,2),'s','MarkerFaceColor',[0.6706 0.7882 0.7176],'MarkerEdgeColor','k','MarkerSize',10)
title('Yearly')

% SE-Deg PW
yearlypCO2_orig = accumarray(year(Deg3.DateTime),Deg3.PorewaterCO2,[],@(x)mean(x,'omitnan'),NaN);

t = zeros(20,2);
count = 1;

for i = 1:10:1000
    y = Deg3((1:i:end),:);
    yearlypCO2 = accumarray(year(y.DateTime),y.PorewaterCO2,[],@(x)mean(x,'omitnan'),NaN);
    if length(yearlypCO2(~isnan(yearlypCO2)))==length(yearlypCO2_orig(~isnan(yearlypCO2_orig)))
        R = corrcoef(yearlypCO2(~isnan(yearlypCO2)&~isnan(yearlypCO2_orig)),yearlypCO2_orig(~isnan(yearlypCO2)&~isnan(yearlypCO2_orig)));
        t(count,1) = 365/i;
        t(count,2) = R(2);
        count = count + 1;
    end
end

z = [z;t];
plot(t(:,1),t(:,2),'o','MarkerFaceColor',[0.7373 0.8196 0.2118],'MarkerEdgeColor',[0.7373 0.8196 0.2118],'MarkerSize',10)

x0 = 0.001;
fitfun = fittype( @(b,x) (x)./(b+x) );
[fitted_curve,~] = fit(z(:,1),z(:,2),fitfun,'StartPoint',x0);
x1 = 0:.1:365;
plot(x1,fitted_curve(x1),'-k','LineWidth',2)

xlim([-2 100])
ylim([0 1.08])
hold off
box on
set(gca,'FontSize',20)
xlabel('Measurement frequency (year^{-1})')
ylabel('Corr. coeff.')
%% Monthly 
nexttile(1)

% US-Myb
monthlypCO2_orig = accumarray(month(USMyb2.DateTime),USMyb2.pCO2,[],@(x)mean(x,'omitnan'),NaN);
t = zeros(9,2);
count = 1;
hold on 

for i = 1:10:1000
    y = USMyb2((1:i:end),:);
    monthlypCO2 = accumarray(month(y.DateTime),y.pCO2,[],@(x)mean(x,'omitnan'),NaN);
    if length(monthlypCO2(~isnan(monthlypCO2)))==12
        R = corrcoef(monthlypCO2(~isnan(monthlypCO2)),monthlypCO2_orig(~isnan(monthlypCO2)));
        t(count,1) = 30/i;
        t(count,2) = R(2);
        count = count + 1;
    end
end

z = t;

plot(t(:,1),t(:,2),'s','MarkerFaceColor',[0.6706 0.7882 0.7176],'MarkerEdgeColor','k','MarkerSize',10)
title('Monthly')

% SE-Deg PW
monthlypCO2_orig = accumarray(month(Deg3.DateTime),Deg3.PorewaterCO2,[],@(x)mean(x,'omitnan'),NaN);
t = zeros(7,2);
count = 1;

for i = 1:10:1000
    y = Deg3((1:i:end),:);
    monthlypCO2 = accumarray(month(y.DateTime),y.PorewaterCO2,[],@(x)mean(x,'omitnan'),NaN);
    if length(monthlypCO2(~isnan(monthlypCO2)))==length(monthlypCO2_orig)
        R = corrcoef(monthlypCO2(~isnan(monthlypCO2)),monthlypCO2_orig(~isnan(monthlypCO2)));
        t(count,1) = 30/i;% on average there are ~30 days per month
        t(count,2) = R(2);
        count = count + 1;
    end
end

z = [z;t];
plot(t(:,1),t(:,2),'o','MarkerFaceColor',[0.7373 0.8196 0.2118],'MarkerEdgeColor',[0.7373 0.8196 0.2118],'MarkerSize',10)

% US-ALQ
monthlypCO2_orig = accumarray(month(USALQ.TIMESTAMP_END),USALQ.pCO2,[],@(x)mean(x,'omitnan'),NaN);
t = zeros(7,2);
count = 1;

for i = 1:10:1000
    y = USALQ((1:i:end),:);
    monthlypCO2 = accumarray(month(y.TIMESTAMP_END),y.pCO2,[],@(x)mean(x,'omitnan'),NaN);
    if length(monthlypCO2(~isnan(monthlypCO2)))==length(monthlypCO2_orig(~isnan(monthlypCO2_orig)))
        R = corrcoef(monthlypCO2(~isnan(monthlypCO2)),monthlypCO2_orig(~isnan(monthlypCO2_orig)));
        t(count,1) = 30/i;
        t(count,2) = R(2);
        count = count + 1;
    end
end

plot(t(:,1),t(:,2),'^','MarkerFaceColor',[0.3333 0.4902 0.0235],'MarkerEdgeColor','k','MarkerSize',10)
z = [z;t];

% US-Los
monthlypCO2_orig = accumarray(month(USLos.TIMESTAMP_END),USLos.pCO2,[],@(x)mean(x,'omitnan'),NaN);
count = 1;
t = zeros(4,2);

for i = 1:10:1000
    y = USLos((1:i:end),:);
    monthlypCO2 = accumarray(month(y.TIMESTAMP_END),y.pCO2,[],@(x)mean(x,'omitnan'),NaN);
    if length(monthlypCO2(~isnan(monthlypCO2)))==length(monthlypCO2_orig(~isnan(monthlypCO2_orig)))
        R = corrcoef(monthlypCO2(~isnan(monthlypCO2)),monthlypCO2_orig(~isnan(monthlypCO2)&~isnan(monthlypCO2_orig)));
        t(count,1) = 30/i;
        t(count,2) = R(2);
        count = count + 1;
    end
end

plot(t(:,1),t(:,2),'o','MarkerFaceColor',[0.6627 0.9882 0.0118],'MarkerEdgeColor','k','MarkerSize',10)
z = [z;t];

x0 = .0001;
fitfun = fittype( @(b,x) (x)./(b+x) );
[fitted_curve,~] = fit(z(:,1),z(:,2),fitfun,'StartPoint',x0);
x1 = 0:.1:16;
plot(x1,fitted_curve(x1),'-k','LineWidth',2)

xlim([-0.5 16])
ylim([0 1.02])
hold off
box on
set(gca,'FontSize',20)
xlabel('Measurement frequency (month^{-1})')
ylabel('Corr. coeff.')
legend({'US-Myb (marsh)','SE-Deg PW (fen)','US-ALQ PW (fen)','US-Los PW (fen)','fit'})

%% Plot daily averages
USLos.pCO2(USLos.pCO2<413)=NaN;
USALQ.pCO2(USALQ.pCO2<415)=NaN;
USMyb2.pCO2(USMyb2.pCO2<405.0)=NaN;
USMyb2 = USMyb2(145:end,:);
USMyb2 = USMyb2(year(USMyb2.DateTime)==2020,:);

nexttile(2)

% US-Myb
t = zeros(21,2);
count = 1;
hold on 

for i = 1:10:1000
    USMybREddyProcout = USMybREddyProcout(year(USMybREddyProcout.DateTime)==2020,:);
    y = USMybREddyProcout((1:i:end),:);
    dailypCO2 = accumarray(day(y.DateTime,'dayofyear'),y.pCO2,[],@(x)mean(x,'omitnan'),NaN);
    if length(dailypCO2)==length(USMyb2.pCO2)
    R = corrcoef(dailypCO2(~isnan(dailypCO2)&~isnan(USMyb2.pCO2)),USMyb2.pCO2(~isnan(dailypCO2)&~isnan(USMyb2.pCO2)));
    t(count,1) = 48/i;
    t(count,2) = R(2);
    count = count + 1;
    end
end

z = t;

plot(t(:,1),t(:,2),'s','MarkerFaceColor',[0.6706 0.7882 0.7176],'MarkerEdgeColor','k','MarkerSize',10)
title('Daily')

% SE-Deg PW
Deg3 = Deg3(year(Deg3.DateTime)==2015,:);
SEDegREddyProcOut2 = SEDegREddyProcOut2(year(SEDegREddyProcOut2.DateTime)==2015,:);

t = zeros(10,2);
count = 1;

for i = 1:10:1000
    y = SEDegREddyProcOut2((1:i:end),:);
    dailypCO2 = accumarray(day(y.DateTime,'dayofyear'),y.pCO2,[],@(x)mean(x,'omitnan'),NaN);
    if length(dailypCO2)==length(Deg3.PorewaterCO2)
    R = corrcoef(dailypCO2(~isnan(dailypCO2)&~isnan(Deg3.PorewaterCO2)),Deg3.PorewaterCO2(~isnan(dailypCO2)&~isnan(Deg3.PorewaterCO2)));
    t(count,1) = 24/i;
    t(count,2) = R(2);
    count = count + 1;
    end
end

z = [z;t];
plot(t(:,1),t(:,2),'o','MarkerFaceColor',[0.7373 0.8196 0.2118],'MarkerEdgeColor',[0.7373 0.8196 0.2118],'MarkerSize',10)

% US-ALQ
t = zeros(20,2);
count = 1;

for i = 1:10:1000
    y = USALQ2((1:i:end),:);
    dailypCO2 = accumarray(day(y.TIMESTAMP_END,'dayofyear'),y.pCO2,[],@(x)mean(x,'omitnan'),NaN);
    if length(dailypCO2)==length(USALQ.pCO2)
        R = corrcoef(dailypCO2(~isnan(dailypCO2)&~isnan(USALQ.pCO2)),USALQ.pCO2(~isnan(dailypCO2)&~isnan(USALQ.pCO2)));
        t(count,1) = 48/i;
        t(count,2) = R(2);
        count = count + 1;
    end
end

plot(t(:,1),t(:,2),'^','MarkerFaceColor',[0.3333 0.4902 0.0235],'MarkerEdgeColor','k','MarkerSize',10)
z = [z;t];

% US-Los
count = 1;
t = zeros(5,2);

for i = 1:10:1000
    y = USLos2((1:i:end),:);
    dailypCO2 = accumarray(day(y.TIMESTAMP_END,'dayofyear'),y.pCO2,[],@(x)mean(x,'omitnan'),0);
    dailypCO2 = dailypCO2(dailypCO2~=0);
    if length(dailypCO2)==length(USLos.pCO2)
        R = corrcoef(dailypCO2(~isnan(dailypCO2)&~isnan(USLos.pCO2)),USLos.pCO2(~isnan(dailypCO2)&~isnan(USLos.pCO2)));
        t(count,1) = 48/i;
        t(count,2) = R(2);
        count = count + 1;
    end
end

plot(t(:,1),t(:,2),'o','MarkerFaceColor',[0.6627 0.9882 0.0118],'MarkerEdgeColor','k','MarkerSize',10)
z = [z;t];

x0 = .0001;
fitfun = fittype( @(b,x) (x)./(b+x) );
[fitted_curve,~] = fit(z(:,1),z(:,2),fitfun,'StartPoint',x0);
x1 = 0:0.001:3;
plot(x1,fitted_curve(x1),'-k','LineWidth',2)

xlim([-0.02 2.02])
ylim([0.8 1.01])
hold off
box on
set(gca,'FontSize',20)
xlabel('Measurement frequency (day^{-1})')

f.TileSpacing = 'compact';
f.Padding = 'compact';
