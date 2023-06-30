%% Loading SubDaily Scale Data
clc;clear;close all;
USLos = readtable("E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls",'FileType','spreadsheet','Sheet','US-Los');
USLos_Daily = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','US-Los_DT');

USALQ2 = readtable("E:\SanDiskSecureAccess\Residencepaper\dielpCO2data.xls",'FileType','spreadsheet','Sheet','US-ALQ_DT');
load USALQSWPCO2.mat
USALQ_PWDaily = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','US-ALQ_DT');
run ALQstreampCO2_2023.m
load USALQStreamDaily.mat

USLos = USLos((ismember(day(USLos.TIMESTAMP_END,'dayofyear'),229:366)&ismember(year(USLos.TIMESTAMP_END),2020)),:);
USLos.pCO2(USLos.pCO2<413)=NaN;
USALQ2.pCO2(USALQ2.pCO2<415)=NaN;
load USEDN.mat
load USMyb.mat;USMyb2 = USMyb2(145:end,:);USMyb2.pCO2(USMyb2.pCO2<405.0)=NaN;
load SEDegPorewater.mat
load SEDegREddyProcOut2.mat
USEDN_Daily = readtable('E:\SanDiskSecureAccess\Residencepaper\pCO2data.xls','Sheet','US-EDN_DT');
%% Calculate correlation between PCO2 and GPP/Reco at each site

% US-EDN
[R,P] = corrcoef(USEDN.CO2(~isnan(USEDN.CO2)&~isnan(USEDN.GPP_DT_uStar)),USEDN.GPP_DT_uStar(~isnan(USEDN.CO2)&~isnan(USEDN.GPP_DT_uStar)));
disp([R P])
[R,P] = corrcoef(USEDN.CO2(~isnan(USEDN.CO2)&~isnan(USEDN.Reco_DT_uStar)),USEDN.Reco_DT_uStar(~isnan(USEDN.CO2)&~isnan(USEDN.Reco_DT_uStar)));
disp([R P])

% US-ALQ porewater 
[R,P] = corrcoef(USALQ_PWDaily.pCO2(~isnan(USALQ_PWDaily.pCO2)&~isnan(USALQ_PWDaily.GPP_F)),USALQ_PWDaily.GPP_F(~isnan(USALQ_PWDaily.pCO2)&~isnan(USALQ_PWDaily.GPP_F)));
disp([R P])
[R,P] = corrcoef(USALQ_PWDaily.pCO2(~isnan(USALQ_PWDaily.pCO2)&~isnan(USALQ_PWDaily.Reco_F)),USALQ_PWDaily.Reco_F(~isnan(USALQ_PWDaily.pCO2)&~isnan(USALQ_PWDaily.Reco_F)));
disp([R P])
% SE-Deg PW
[R,P] = corrcoef(SEDegREddyProcOut2.pCO2(~isnan(SEDegREddyProcOut2.pCO2)&~isnan(SEDegREddyProcOut2.GPP_DT_uStar)),SEDegREddyProcOut2.GPP_DT_uStar(~isnan(SEDegREddyProcOut2.pCO2)&~isnan(SEDegREddyProcOut2.GPP_DT_uStar)));
disp([R P])
[R,P] = corrcoef(SEDegREddyProcOut2.pCO2(~isnan(SEDegREddyProcOut2.pCO2)&~isnan(SEDegREddyProcOut2.Reco_DT_uStar)),SEDegREddyProcOut2.Reco_DT_uStar(~isnan(SEDegREddyProcOut2.pCO2)&~isnan(SEDegREddyProcOut2.Reco_DT_uStar)));
disp([R P])
%% Plotting
[R,P] = corrcoef(USLos.pCO2(~isnan(USLos.pCO2)&~isnan(USLos.TA_PI_F_1_1_1)),USLos.TA_PI_F_1_1_1(~isnan(USLos.pCO2)&~isnan(USLos.TA_PI_F_1_1_1)));
disp([R P])
    
% US-Los
% Calculate sensitivity parameter
q = zeros(1,51);
count = 1;

for b = 0:0.001:0.05
    pCO2norm = USLos_Daily.pCO2.*(exp(b*(mean(USLos_Daily.Tair_f,'omitnan')-USLos_Daily.Tair_f)));
    x2 = corrcoef(pCO2norm(~isnan(pCO2norm)&~isnan(USLos_Daily.Tair_f)),USLos_Daily.Tair_f(~isnan(pCO2norm)&~isnan(USLos_Daily.Tair_f)));
    q(count) = x2(2);
    count = count+1;
end

if any(q>0)
    ind = find(q==min(q(q>=0)));
else
    ind = find(q==max(q(q<=0)));
end

b = 0:0.001:0.05;

% Remove influence of temperature with new calibrated sensitivity parameter
USLos.pCO2 = USLos.pCO2.*(exp(b(ind)*(mean(USLos.TA_PI_F_1_1_1,'omitnan')-USLos.TA_PI_F_1_1_1)));
[R,P] = corrcoef(USLos.pCO2(~isnan(USLos.pCO2)&~isnan(USLos.TA_PI_F_1_1_1)),USLos.TA_PI_F_1_1_1(~isnan(USLos.pCO2)&~isnan(USLos.TA_PI_F_1_1_1)));
disp([R P])

%% US-Myb
[R,P] = corrcoef(USMybREddyProcout.pCO2(~isnan(USMybREddyProcout.pCO2)&~isnan(USMybREddyProcout.Twater)),USMybREddyProcout.Twater(~isnan(USMybREddyProcout.pCO2)&~isnan(USMybREddyProcout.Twater)));
disp([R P])

% Calculate sensitivity parameter
q = zeros(1,51);
count = 1;

for b = 0:0.001:0.05
    pCO2norm = USMyb2.pCO2.*(exp(b*(mean(USMyb2.Twater,'omitnan')-USMyb2.Twater)));
    x2 = corrcoef(pCO2norm(~isnan(pCO2norm)&~isnan(USMyb2.Twater)),USMyb2.Twater(~isnan(pCO2norm)&~isnan(USMyb2.Twater)));
    q(count) = x2(2);
    count = count+1;
end

if any(q>0)
    ind = find(q==min(q(q>=0)));
else
    ind = find(q==max(q(q<=0)));
end

b = 0:0.001:0.05;
% Remove influence of temperature with new calibrated sensitivity parameter
USMybREddyProcout.pCO2 = USMybREddyProcout.pCO2.*(exp(b(ind)*(mean(USMybREddyProcout.Twater,'omitnan')-USMybREddyProcout.Twater)));
[R,P] = corrcoef(USMybREddyProcout.pCO2(~isnan(USMybREddyProcout.pCO2)&~isnan(USMybREddyProcout.Twater)),USMybREddyProcout.Twater(~isnan(USMybREddyProcout.pCO2)&~isnan(USMybREddyProcout.Twater)));
disp([R P])
%% SE-Deg PW
[R,P] = corrcoef(SEDegREddyProcOut2.pCO2(~isnan(SEDegREddyProcOut2.pCO2)&~isnan(SEDegREddyProcOut2.Twater)),SEDegREddyProcOut2.Twater(~isnan(SEDegREddyProcOut2.pCO2)&~isnan(SEDegREddyProcOut2.Twater)));
disp([R P])

% Calculate sensitivity parameter
q = zeros(1,51);
count = 1;

for b = 0:0.001:0.05
    pCO2norm = Deg3.PorewaterCO2.*(exp(b*(mean(Deg3.Twater,'omitnan')-Deg3.Twater)));
    x2 = corrcoef(pCO2norm(~isnan(pCO2norm)&~isnan(Deg3.Twater)),Deg3.Twater(~isnan(pCO2norm)&~isnan(Deg3.Twater)));
    q(count) = x2(2);
    count = count+1;
end

if any(q>0)
    ind = find(q==min(q(q>=0)));
else
    ind = find(q==max(q(q<=0)));
end

b = 0:0.001:0.05;

% Remove influence of temperature with new calibrated sensitivity parameter
SEDegREddyProcOut2.pCO2 = SEDegREddyProcOut2.pCO2.*(exp(b(ind)*(mean(SEDegREddyProcOut2.Twater,'omitnan')-SEDegREddyProcOut2.Twater)));
[R,P] = corrcoef(SEDegREddyProcOut2.pCO2(~isnan(SEDegREddyProcOut2.pCO2)&~isnan(SEDegREddyProcOut2.Twater)),SEDegREddyProcOut2.Twater(~isnan(SEDegREddyProcOut2.pCO2)&~isnan(SEDegREddyProcOut2.Twater)));
disp([R P])
%% USALQ PW
[R,P] = corrcoef(USALQ2.pCO2(~isnan(USALQ2.pCO2)&~isnan(USALQ2.TA_F)),USALQ2.TA_F(~isnan(USALQ2.pCO2)&~isnan(USALQ2.TA_F)));
disp([R P])
% Calculate sensitivity parameter
q = zeros(1,51);
count = 1;

for b = 0:0.001:0.05
    pCO2norm = USALQ_PWDaily.pCO2.*(exp(b*(mean(USALQ_PWDaily.TA_F,'omitnan')-USALQ_PWDaily.TA_F)));
    x2 = corrcoef(pCO2norm(~isnan(pCO2norm)&~isnan(USALQ_PWDaily.TA_F)),USALQ_PWDaily.TA_F(~isnan(pCO2norm)&~isnan(USALQ_PWDaily.TA_F)));
    q(count) = x2(2);
    count = count+1;
end

if any(q>0)
    ind = find(q==min(q(q>=0)));
else
    ind = find(q==max(q(q<=0)));
end

b = 0:0.001:0.05;

% Remove influence of temperature with new calibrated sensitivity parameter
USALQ2.pCO2 = USALQ2.pCO2.*(exp(b(ind)*(mean(USALQ2.TA_F,'omitnan')-USALQ2.TA_F)));
[R,P] = corrcoef(USALQ2.pCO2(~isnan(USALQ2.pCO2)&~isnan(USALQ2.TA_F)),USALQ2.TA_F(~isnan(USALQ2.pCO2)&~isnan(USALQ2.TA_F)));
disp([R P])
%% US-ALQ SW
[R,P] = corrcoef(USALQ.SWpCO2(~isnan(USALQ.SWpCO2)&~isnan(USALQ.TA_F)),USALQ.TA_F(~isnan(USALQ.SWpCO2)&~isnan(USALQ.TA_F)));
disp([R P])
% Calculate sensitivity parameter
q = zeros(1,51);
count = 1;

for b = 0:0.001:0.05 
    pCO2norm = USALQStreamDaily.pCO2_ppm.*(exp(b*(mean(USALQStreamDaily.Water_Temp_C,'omitnan')-USALQStreamDaily.Water_Temp_C)));
    x2 = corrcoef(pCO2norm(~isnan(pCO2norm)&~isnan(USALQStreamDaily.Water_Temp_C)),USALQStreamDaily.Water_Temp_C(~isnan(pCO2norm)&~isnan(USALQStreamDaily.Water_Temp_C)));
    q(count) = x2(2);
    count = count+1;
end

if any(q>0)
    ind = find(q==min(q(q>=0)));
else
    ind = find(q==max(q(q<=0)));
end

b = 0:0.001:0.05;

% Remove influence of temperature with new calibrated sensitivity parameter
USALQ.SWpCO2 = USALQ.SWpCO2.*(exp(b(ind)*(mean(USALQ.TA_F,'omitnan')-USALQ.TA_F)));
[R,P] = corrcoef(USALQ.SWpCO2(~isnan(USALQ.SWpCO2)&~isnan(USALQ.TA_F)),USALQ.TA_F(~isnan(USALQ.SWpCO2)&~isnan(USALQ.TA_F)));
disp([R P])
%% US-EDN
[R,P] = corrcoef(USEDN.CO2(~isnan(USEDN.CO2)&~isnan(USEDN.Wtemp)),USEDN.Wtemp(~isnan(USEDN.CO2)&~isnan(USEDN.Wtemp)));
disp([R P])

% Calculate sensitivity parameter
q = zeros(1,51);
count = 1;

for b = 0:0.001:0.05
    pCO2norm = USEDN_Daily.CO2.*(exp(b*(mean(USEDN_Daily.Wtemp,'omitnan')-USEDN_Daily.Wtemp)));
    x2 = corrcoef(pCO2norm(~isnan(pCO2norm)&~isnan(USEDN_Daily.Wtemp)),USEDN_Daily.Wtemp(~isnan(pCO2norm)&~isnan(USEDN_Daily.Wtemp)));
    q(count) = x2(2);
    count = count+1;
end

if any(q>0)
    ind = find(q==min(q(q>=0)));
else
    ind = find(q==max(q(q<=0)));
end

b = 0:0.001:0.05;

USEDN.CO2 = USEDN.CO2.*(exp(b(ind)*(mean(USEDN.Wtemp,'omitnan')-USEDN.Wtemp)));
[R,P] = corrcoef(USEDN.CO2(~isnan(USEDN.CO2)&~isnan(USEDN.Wtemp)),USEDN.Wtemp(~isnan(USEDN.CO2)&~isnan(USEDN.Wtemp)));
disp([R P])
%% Plotting US-ALQ PW
dailypCO2_min = accumarray(day(USALQ2.TIMESTAMP_END,'dayofyear'),USALQ2.pCO2,[],@min,NaN);
dailypCO2_max = accumarray(day(USALQ2.TIMESTAMP_END,'dayofyear'),USALQ2.pCO2,[],@max,NaN);
dailypCO2_bioALQ = dailypCO2_max-dailypCO2_min;
dailyGPPALQ = accumarray(day(USALQ2.TIMESTAMP_END,'dayofyear'),USALQ2.RECO_F,[],@(x)mean(x,'omitnan'),NaN);

figure()
plot(dailypCO2_bioALQ,dailyGPPALQ,'.')

um = USALQ2.pCO2(~isnan(USALQ2.pCO2)&~isnan(USALQ2.RECO_F));
un = USALQ2.RECO_F(~isnan(USALQ2.pCO2)&~isnan(USALQ2.RECO_F));
fs = 48/86400;
figure()
mscohere(um,un,hamming(70),[],[],fs);

h = gcf;
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children');
xdata = get(dataObjs, 'XData'); 
ydata = get(dataObjs, 'YData');

%% Repeat with US-Los
dailypCO2_min = accumarray(day(USLos.TIMESTAMP_END,'dayofyear'),USLos.pCO2,[],@min,NaN);
dailypCO2_max = accumarray(day(USLos.TIMESTAMP_END,'dayofyear'),USLos.pCO2,[],@max,NaN);
dailypCO2_bioLos = dailypCO2_max-dailypCO2_min;
dailyGPPLos = accumarray(day(USLos.TIMESTAMP_END,'dayofyear'),USLos.RECO_PI_F,[],@(x)mean(x,'omitnan'),NaN);

figure()
plot(dailypCO2_bioLos,dailyGPPLos,'.')

um = USLos.pCO2(~isnan(USLos.pCO2)&~isnan(USLos.RECO_PI_F));
un = USLos.RECO_PI_F(~isnan(USLos.pCO2)&~isnan(USLos.RECO_PI_F));
fs = 48/86400;

figure()
mscohere(um,un,hamming(70),[],[],fs);

h = gcf;
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children');
xdata2 = get(dataObjs, 'XData'); 
ydata2 = get(dataObjs, 'YData');

%% Repeat with US-Myb
dailypCO2_min = accumarray(day(USMybREddyProcout.DateTime,'dayofyear'),USMybREddyProcout.pCO2,[],@min,NaN);
dailypCO2_max = accumarray(day(USMybREddyProcout.DateTime,'dayofyear'),USMybREddyProcout.pCO2,[],@max,NaN);
dailypCO2_bioMyb = dailypCO2_max-dailypCO2_min;
dailyGPPMyb = accumarray(day(USMybREddyProcout.DateTime,'dayofyear'),USMybREddyProcout.Reco_DT_uStar,[],@(x)mean(x,'omitnan'),NaN);

figure()
plot(dailypCO2_bioMyb,dailyGPPMyb,'.')

um = USMybREddyProcout.pCO2(~isnan(USMybREddyProcout.pCO2)&~isnan(USMybREddyProcout.Reco_DT_uStar));
un = USMybREddyProcout.Reco_DT_uStar(~isnan(USMybREddyProcout.pCO2)&~isnan(USMybREddyProcout.Reco_DT_uStar));
fs = 48/86400;

figure()
mscohere(um,un,hamming(70),[],[],fs);

h = gcf;
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children');
xdata3 = get(dataObjs, 'XData'); 
ydata3 = get(dataObjs, 'YData');

%% Repeat with SE-Deg
dailypCO2_min = accumarray(day(SEDegREddyProcOut2.DateTime,'dayofyear'),SEDegREddyProcOut2.pCO2,[],@min,NaN);
dailypCO2_max = accumarray(day(SEDegREddyProcOut2.DateTime,'dayofyear'),SEDegREddyProcOut2.pCO2,[],@max,NaN);
dailypCO2_bioDeg = dailypCO2_max-dailypCO2_min;
dailyGPPDeg = accumarray(day(SEDegREddyProcOut2.DateTime,'dayofyear'),SEDegREddyProcOut2.Reco_DT_uStar,[],@(x)mean(x,'omitnan'),NaN);

figure()
plot(dailypCO2_bioDeg,dailyGPPDeg,'.')

um = SEDegREddyProcOut2.pCO2(~isnan(SEDegREddyProcOut2.pCO2)&~isnan(SEDegREddyProcOut2.Reco_DT_uStar));
un = SEDegREddyProcOut2.Reco_DT_uStar(~isnan(SEDegREddyProcOut2.pCO2)&~isnan(SEDegREddyProcOut2.Reco_DT_uStar));
fs = 48/86400;

figure()
mscohere(um,un,hamming(70),[],[],fs);

h = gcf;
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children');
xdata4 = get(dataObjs, 'XData'); 
ydata4 = get(dataObjs, 'YData');
%% Repeat with US-EDN
dailypCO2_min = accumarray(day(USEDN.DateTime,'dayofyear'),USEDN.CO2,[],@min,NaN);
dailypCO2_max = accumarray(day(USEDN.DateTime,'dayofyear'),USEDN.CO2,[],@max,NaN);
dailypCO2_bioEDN = dailypCO2_max-dailypCO2_min;
dailyGPPEDN = accumarray(day(USEDN.DateTime,'dayofyear'),USEDN.Reco_DT_uStar,[],@(x)mean(x,'omitnan'),NaN);

figure()
plot(dailypCO2_bioEDN,dailyGPPEDN,'.')

um = USEDN.CO2(~isnan(USEDN.CO2)&~isnan(USEDN.Reco_DT_uStar));
un = USEDN.Reco_DT_uStar(~isnan(USEDN.CO2)&~isnan(USEDN.Reco_DT_uStar));
fs = 48/86400;

figure()
mscohere(um,un,hamming(70),[],[],fs);

h = gcf;
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children');
xdata5 = get(dataObjs, 'XData'); 
ydata5 = get(dataObjs, 'YData');
%% Repeat with US-ALQ Stream

dailypCO2_min = accumarray(day(USALQ.TIMESTAMP_END,'dayofyear'),USALQ.SWpCO2,[],@min,NaN);
dailypCO2_max = accumarray(day(USALQ.TIMESTAMP_END,'dayofyear'),USALQ.SWpCO2,[],@max,NaN);
dailypCO2_bioALQSW = dailypCO2_max-dailypCO2_min;
dailyGPPALQSW = accumarray(day(USALQ.TIMESTAMP_END,'dayofyear'),USALQ.RECO_F,[],@(x)mean(x,'omitnan'),NaN);

um = USALQ.SWpCO2(~isnan(USALQ.SWpCO2)&~isnan(USALQ.RECO_F));
un = USALQ.GPP_F(~isnan(USALQ.SWpCO2)&~isnan(USALQ.RECO_F));
fs = 48/86400;

figure()
mscohere(um,un,hamming(70),[],[],fs);

h = gcf;
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children');
xdata6 = get(dataObjs, 'XData'); 
ydata6 = get(dataObjs, 'YData');

% close open figure windows
close all
%% Plot all together
dailypCO2_bioDeg(129) = NaN;

binEdges = 0:1000:2e4;
groupPCO2BIO = discretize(dailypCO2_bioMyb,binEdges);
[meanGPP,BG] = groupsummary(dailyGPPMyb,groupPCO2BIO,'mean');

figure()
z = tiledlayout(1,3);
nexttile(2)
hold on
plot(binEdges(BG),meanGPP,'d','MarkerFaceColor',[0.4627 0.6392 0.5804],'MarkerEdgeColor','k','MarkerSize',12)
coefficients = polyfit(binEdges(BG),meanGPP,1);
[R,P] = corrcoef(binEdges(BG),meanGPP);
disp([R P])

xFit = binEdges;
yFit = polyval(coefficients,xFit);
plot(xFit, yFit,'-','Color',[0.4627 0.6392 0.5804],'LineWidth',4); % Plot fitted line.
hold off 
legend('US-Myb SW (Marsh)','','Location','northwest')
xlim([0 12000])
set(gca,'FontSize',22)
ylim([0 12])
box on
xlabel('\DeltaP_{CO2 at Tmean}')
yticklabels({''})

nexttile(1)
hold on
binEdges = 0:50:3000;
groupPCO2BIO = discretize(dailypCO2_bioEDN,binEdges);
[meanGPPEDN,BGEDN] = groupsummary(dailyGPPEDN,groupPCO2BIO,'mean');
plot(binEdges(BGEDN(~isnan(BGEDN))),meanGPPEDN(~isnan(BGEDN)),'o','MarkerFaceColor',[0.7882 0.6000 1.0000],'MarkerEdgeColor','k','MarkerSize',12)

coefficients = polyfit(binEdges(BGEDN(~isnan(BGEDN))),meanGPPEDN(~isnan(BGEDN)),1);
xFit = binEdges;
yFit = polyval(coefficients,xFit);
plot(xFit, yFit,'-','Color',[0.7882 0.6000 1.0000],'LineWidth',4); % Plot fitted line.

groupPCO2BIO = discretize(dailypCO2_bioALQSW,binEdges);
[meanGPPALQSW,BGALQSW] = groupsummary(dailyGPPALQSW,groupPCO2BIO,'mean');
plot(binEdges(BGALQSW(~isnan(BGALQSW))),meanGPPALQSW(~isnan(BGALQSW)),'v','MarkerFaceColor',[0.3333 0.4902 0.0235],'MarkerEdgeColor','k','MarkerSize',12)

coefficients = polyfit(binEdges(BGALQSW(~isnan(BGALQSW))),meanGPPALQSW(~isnan(BGALQSW)),1);
xFit = binEdges;
yFit = polyval(coefficients,xFit);
plot(xFit, yFit,'-','Color',[0.3333 0.4902 0.0235],'LineWidth',4); % Plot fitted line.

hold off
legend('US-EDN SW (Tidal)','','US-ALQ SW (Fen)','','Location','northwest')
ylim([0 12])
ylabel('R_{eco} (\mumol CO_2 m^-^2s^-^1)')
set(gca,'FontSize',22)
box on

binEdges = 0:5000:5e4;
groupPCO2BIO = discretize(dailypCO2_bioALQ,binEdges);
[meanGPPALQ,BGALQ] = groupsummary(dailyGPPALQ,groupPCO2BIO,'mean');

coefficients = polyfit(binEdges(BGALQ(~isnan(BGALQ))),meanGPPALQ(~isnan(BGALQ)),1);
xFitALQ = 0:5000:10e4;
yFitALQ = polyval(coefficients,xFitALQ);
[R,P] = corrcoef(binEdges(BGALQ(~isnan(BGALQ))),meanGPPALQ(BGALQ(~isnan(BGALQ))));
disp([R P])

groupPCO2BIO = discretize(dailypCO2_bioLos,binEdges);
[meanGPPLos,BGLos] = groupsummary(dailyGPPLos,groupPCO2BIO,'mean');
[R,P] = corrcoef(binEdges(BGLos(~isnan(BGLos))),meanGPPLos(BGLos(~isnan(BGLos))));
disp([R P])

coefficients = polyfit(binEdges(BGLos(~isnan(BGLos))),meanGPPLos(~isnan(BGLos)),1);

xFitLos = 0:5000:10e4;
yFitLos = polyval(coefficients,xFitLos);

ax = nexttile(3);
hold on
plot(binEdges(BGALQ(~isnan(BGALQ))),meanGPPALQ(~isnan(BGALQ)),'o','MarkerFaceColor',[0.3333 0.4902 0.0235],'MarkerEdgeColor','k','MarkerSize',12)
plot(xFitALQ, yFitALQ,'-','Color',[0.3333 0.4902 0.0235],'LineWidth',4); % Plot fitted line.

plot(binEdges(BGLos(~isnan(BGLos))),meanGPPLos(~isnan(BGLos)),'s','MarkerFaceColor',[0.1843 0.7608 0.1843],'MarkerEdgeColor','k','MarkerSize',12)
plot(xFitLos, yFitLos,'-','Color',[0.1843 0.7608 0.1843],'LineWidth',4); % Plot fitted line.

groupPCO2BIO = discretize(dailypCO2_bioDeg,binEdges);
[meanGPPDeg,BGDeg] = groupsummary(dailyGPPDeg,groupPCO2BIO,'mean');
plot(binEdges(BGDeg(~isnan(BGDeg))),meanGPPDeg(~isnan(BGDeg)),'^','MarkerFaceColor',[0.7373 0.8196 0.2118],'MarkerEdgeColor','k','MarkerSize',12)
coefficients = polyfit(binEdges(BGDeg(~isnan(BGDeg))),meanGPPDeg(~isnan(BGDeg)),1);
[R,P] = corrcoef(binEdges(BGDeg(~isnan(BGDeg))),meanGPPDeg(~isnan(BGDeg)));
disp([R P])

xFit = binEdges;
yFit = polyval(coefficients,xFit);
plot(xFit, yFit,'-','Color',[0.7373 0.8196 0.2118],'LineWidth',4); % Plot fitted line.

hold off
ylim([0 12])
xlim([0 40000])
ylabel('R_{eco} (\mumol CO_2 m^-^2s^-^1)')
set(gca,'FontSize',22)
legend('US-ALQ PW (fen)','','US-Los PW (fen)','','SE-Deg PW (fen)','')
box on

z.TileSpacing = 'compact';
z.Padding = 'compact';
ax.YAxisLocation = 'right';
