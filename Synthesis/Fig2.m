clear; clc;
T = readtable('E:\SanDiskSecureAccess\pCO2data\SiteLocations.xlsx');
GT = table2geotable(T,'geographic');
latlim = [20 85];lonlim = [-180 30];
ax = worldmap(latlim,lonlim);
load koppenlegend.mat;

sitenums = string(T.Var1);
sitenums(1,1) = "1 & 2 ";
sitenums(2,1) = "";
sitenums(3,1) = "3  ";
sitenums(4,1) = "  ";
sitenums(5,1) = "5  ";
sitenums(6,1) = "6 & 7 ";
sitenums(7,1) = "";
sitenums(8,1) = "8  ";
sitenums(9,1) = "9 & 10  ";
sitenums(10,1) = "";
sitenums(11,1) = "";
sitenums(12,1) = "11 & 12  ";
sitenums(13,1) = "13 & 22  ";
sitenums(14,1) = "14 ";
sitenums(15,1) = "15 & 19 ";
sitenums(16,1) = "16 & 4  ";
sitenums(17,1) = "17  ";
sitenums(18,1) = "18  ";
sitenums(19,1) = "  ";
sitenums(20,1) = "20  ";
sitenums(21,1) = "21  ";
sitenums(22,1) = "  ";

textm(T.Latitude,T.Longitude,sitenums, ...
  "HorizontalAlignment","right",'FontSize',20,'fontweight','bold','Color','k')
mlabel("off")
load coastlines
plotm(coastlat,coastlon,'Color','w')

wetlandtypes = string(T.WetlandType);
fens = find(ismember(wetlandtypes,'fen')|ismember(wetlandtypes,'fen & bog'));
tidal = find(ismember(wetlandtypes,'tidal'));
bog = find(ismember(wetlandtypes,'bog'));
marsh = find(ismember(wetlandtypes,'marsh'));
prairies = find(ismember(wetlandtypes,'prairie pothole or karst'));
alpine = find(ismember(wetlandtypes,'alpine'));

geoshow(ax,GT(fens,:),'Marker','p','MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',2);
geoshow(ax,GT(tidal,:),'Marker','*','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',2);
geoshow(ax,GT(bog,:),'Marker','diamond','MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','w','LineWidth',2);
geoshow(ax,GT(marsh,:),'Marker','square','MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','w','LineWidth',2);
geoshow(ax,GT(prairies,:),'Marker','o','MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','w','LineWidth',2);
geoshow(ax,GT(alpine,:),'Marker','^','MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','w','LineWidth',2);

[A,R] = readgeoraster("C:\Users\lturn\Downloads\Beck_KG_V1\Beck_KG_V1_present_0p5.tif", 'OutputType', 'double');
h = geoshow(A,R,"DisplayType","surface", 'FaceAlpha', 0.8);
uistack(h,'bottom')
setm(gca,'FontSize',15)

% add legend
newStr = split(legend.VarName4,[" ","]","["]);newStr(:,1) = [];
newStr = str2double(newStr)./256;
newStr(:,4) = [];
cmap = newStr;
c = colorbar('Ticks',1:1:30,'TickLabels',legend.Af,'FontSize',15,'Location','southoutside');
