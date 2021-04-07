%% simulation paper - get bimodal region
clear
close all
clc
seed = 'seed2';

region = importdata(['/Volumes/Seagate Backup Plus Drive/ForViraNewResponders/KATP_MovetoAvg_' seed '/ZeroCoupCell.txt']);
region = region+1;

positions = importdata(['/Volumes/Seagate Backup Plus Drive/ForViraNewResponders/KATP_MovetoAvg_' seed '/XYZpos.txt']);
x = positions(:,1);
y = positions(:,2);
z = positions(:,3);

c = zeros(1,size(positions,1));
c(region)=1;
maxc = max(c);
%c = c/maxc;

%%size of circles%%
s = repmat(150,1000,1);

bwr = @(n)interp1([1 2 3], [0 0 1; 1 1 1; 1 0  0], linspace(1, 3, n), 'linear')
%colormap(bwr(64));
cmap = bwr(64);
% cmap = gray(10);
% cmap = cmap(4:end,:);

% cmap = parula(4);
% cmap(4, :) = 1;

cmap = [1 0  0; 1 1 1];

figure()
h = scatter3(x,y,z,s, c);
h.MarkerFaceColor = 'flat';
h.MarkerEdgeColor = 'k';
colormap(flipud(cmap));
colorbar
caxis([0 1])
grid off
set(gca,'XLim',[-.5 1.5],'YLim',[-.5 1.5],'ZLim',[-.5 1.5])


%% 
region2 = importdata('/Users/jdwulet/Google Drive/LabStuff/Projects/FirstResponder/IndividualData/KATPSeed2FirstResponder100.txt');
region2 = region2+1;
%region2 = sort(region2);

 memofnewregion = ismember(region2, region);
 overlap = region2.*memofnewregion;
 getoverlapind = find(overlap > 0);
 justoverlap = overlap(getoverlapind);
  
 notmemofnewregion = ~ismember(region2, region);
 notoverlap = region2.*notmemofnewregion;
newind = find(notoverlap > 0);
 onlynew = notoverlap(newind);
 
  origmem = ~ismember(region, region2);
 overlap2 = region.*origmem;
 oldind = find(overlap2 > 0);
 onlyold = overlap2(oldind);
%region = region;

c = zeros(1,size(positions,1));
c(region)=1;
c(region2)=2;
c(justoverlap)=3;
maxc = max(c);
%c = c/maxc;

%%size of circles%%
s = repmat(150,1000,1);

cmap = [ 1 0  0; 1 .6  .6; .5 .5 .5; 1 1 1];

figure()
h = scatter3(x,y,z,s, c);
h.MarkerFaceColor = 'flat';
h.MarkerEdgeColor = 'k';
colormap(flipud(cmap));
colorbar
caxis([0 3])
grid off
set(gca,'XLim',[-.5 1.5],'YLim',[-.5 1.5],'ZLim',[-.5 1.5])

%% get calcium traces
calciumorig = importdata(['/Volumes/Seagate Backup Plus Drive/ForViraNewResponders/WT_G2to11_Coup120Seed2/calcium.txt']);
calciumnew = importdata(['/Volumes/Seagate Backup Plus Drive/ForViraNewResponders/KATP_MovetoAvg_' seed '/calcium.txt']);

averageorig = mean(calciumorig, 2);
averagenew = mean(calciumnew, 2);
averagenew = averagenew(1:999);

ninecells = [onlyold(1:3); justoverlap(1:3); onlynew(1:3)];
origtc = calciumorig(:,ninecells);
newtc= calciumnew(:,ninecells);
newtc = newtc(1:999, :);

figure
plot(origtc);
legend

figure
plot(newtc);

%% get first, last and average responder
calciumorig = importdata(['/Volumes/Seagate Backup Plus Drive/ForViraNewResponders/WT_G2to11_Coup120Seed1/calcium.txt']);
responderlist = importdata(['/Users/jdwulet/Google Drive/LabStuff/Projects/FirstResponder/IndividualData/WT100Seed1FullListofResponder.txt']);

averageorig = mean(calciumorig, 2);
first = calciumorig(:,responderlist(1:3));
last = calciumorig(:,responderlist(1000-3:1000));

figure
plot(averageorig,'k');
hold on
plot(first, 'r');
hold on
plot(last, 'b');
legend('avg', 'first', 'last');

%%get first responder time scross islet
%% simulation paper %phase lag continuous

c =  importdata('/Users/jdwulet/Google Drive/LabStuff/Projects/FirstResponder/WT100Seed1ResponseTime.txt');

positions = importdata(['/Volumes/Seagate Backup Plus Drive/ForViraNewResponders/WT_G2to11_Coup120Seed1/XYZpos.txt']);
x = positions(:,1);
y = positions(:,2);
z = positions(:,3);

%%size of circles%%
s = repmat(150,1000,1);

bwr = @(n)interp1([1 2 3], [0 0 1; 1 1 1; 1 0  0], linspace(1, 3, n), 'linear');
%colormap(bwr(64));

Fcs = [0.8,0,0]; %RespCell colormap start
Fce = [0,0.8,0.1]; %RespCell colormap end
FstR = linspace(Fcs(1),Fce(1))';
FstG = linspace(Fcs(2),Fce(2))';
FstB = linspace(Fcs(3),Fce(3))';
FirstMap = [FstR,FstG,FstB]; %Generates colormap for RespCells


figure
h = scatter3(x,y,z,s, c);
h.MarkerFaceColor = 'flat';
h.MarkerEdgeColor = 'k';
colormap(FirstMap)
%colormap(flipud(bwr(64)));
%colormap(bwr(64));
colorbar
%caxis([0.16 0.2])
grid off
set(gca,'XLim',[-.5 1.5],'YLim',[-.5 1.5],'ZLim',[-.5 1.5])
title('respnsetime')

