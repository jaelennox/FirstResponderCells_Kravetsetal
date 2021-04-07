clear all
close all
clc

savename = 'checkuncoupled_perc10';
name = {'uncoupled_perc10_v2', 'uncoupled_random_perc10_v2'};
val = 100;

savename = 'checkuncoupled_perc05';
name = {'uncoupled_perc05_v2', 'uncoupled_random_perc05_v2'};
val = 5;
% savename = 'checkuncoupled_perc1';
% name = {'uncoupled_perc1_v2', 'uncoupled_random_perc1_v2'};
% val = 10;
% savename = 'checkuncoupled_perc20';
% name = {'uncoupled_perc20_v2', 'uncoupled_random_perc20_v2'};
% val = 200;
% savename = 'checkuncoupled_perc30';
% name = {'uncoupled_perc30_v2', 'uncoupled_random_perc30_v2'};
% val = 300;
% savename = 'checkuncoupled_perc50';
% name = {'uncoupled_perc50_v2', 'uncoupled_random_perc50_v2'};
% val = 500;


seed = 1:5;
count = 1;

v=1;
mainpath = '/Volumes/Seagate Backup Plus Drive/ForViraNewResponders/';
resultspath = '/Users/jdwulet/Google Drive/LabStuff/Projects/FirstResponder/';

for i = 1:length(name)
    for ii = 1:length(seed)
        origfull = importdata([resultspath 'WT100Seed' num2str(seed(ii)) 'FullListofResponder.txt']);
        
        get100values = origfull(val+1:val+100);
        
        newfirst = importdata([resultspath name{i} 'Seed' num2str(seed(ii)) 'FirstResponder' num2str(val) '.txt']);
        %secondresponder = importdata([name{i} 'Seed' num2str(seed(ii)) 'SecondResponder.txt']);
        
        memoffirst = ismember(newfirst, get100values);
          % memofsecond = ismember(newfirst, origsecond);
        
        Data.name{count, 1} = name{i};
        Data.seed(count,1) = seed(ii);
        Data.percentfirst(count,1) = sum(memoffirst)/length(get100values);
        % Data.percentsecond(count,1) = sum(memofsecond)/length(origfirst);
        count=count+1;
        
    end
end

filename = [resultspath savename datestr(datetime('today'))];
T = struct2table(Data);
writetable(T,[filename '.xlsx']);
%