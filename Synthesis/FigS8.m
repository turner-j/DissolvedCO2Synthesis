%% Clear workspace and command window and load data
clear ;clc;close all; 
addpath 'E:\SanDiskSecureAccess\all MATLAB files\allequash'
cd('E:\SanDiskSecureAccess\all MATLAB files\Residence paper\UsEvM')
load('USEvm_pCO2_clean.mat')
run ImportUSEvMREddyProcOut.m
%% Index for similar dates and plot
t1 = dateshift(SALOGR08.datetime(1),'start','hour');
t2 = dateshift(SALOGR08.datetime(end),'end','hour');

ind = find(ismember(USEvMREddyProcOut.DateTime,t1));
ind2 = find(ismember(USEvMREddyProcOut.DateTime,t2));
USEvMREddyProcOut = USEvMREddyProcOut((ind:ind2),:);

%% Diel pCO2 cycles
hourlypCO2 = accumarray(hour(SALOGR08.datetime)+1,SALOGR08.CO2ppm,[],@nanmean);
hrs = 0:0.5:23.5;

hourlyTwater= accumarray(hour(USEvMREddyProcOut.DateTime)+1,USEvMREddyProcOut.SoilWTemp,[],@nanmean);

figure()
hold on 
ix = find(~ismember(unique(hour(USEvMREddyProcOut.DateTime)+1),9:16));% Nighttime is anything NOT 8am-3pm
scatter(hourlyTwater(ix),hourlypCO2(ix),50,'filled','MarkerFaceColor','k','LineWidth',1.5,'MarkerEdgeColor','k')
scatter(hourlyTwater,hourlypCO2,'LineWidth',1.5,'MarkerEdgeColor','k')
hold off
xlabel('T_{water} ({\circ}C)')
ylabel('CO_2 {\itaq} (ppmv)')
set(gca,'FontSize',15)
ylim([750 3000])
xlim([28 33.5])
legend({'Nighttime','Daytime'})
