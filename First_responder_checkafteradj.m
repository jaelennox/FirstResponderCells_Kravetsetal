clear all
close all
clc

savename = 'checkafteradj';
name = {'Coup', 'GK', 'KATP'};

% savename = 'checkafteradj_XYZ';
% name = {'WT100_changeXYZ', 'WT100_changeXYZ_newsphere'};

seed = 1:5;
count = 1;
val = [100];
v=1;
mainpath = '/Users/jdwulet/Google Drive/LabStuff/Projects/FirstResponder/';
resultspath = '/Users/jdwulet/Google Drive/LabStuff/Projects/FirstResponder/';

for i = 1:length(name)
    for ii = 1:length(seed)
        
        origfirst = importdata([resultspath 'WT100' 'Seed' num2str(seed(ii)) 'FirstResponder' num2str(val(v)) '.txt']);
        origsecond = importdata([resultspath 'WT100' 'Seed' num2str(seed(ii)) 'SecondResponder' num2str(val(v)) '.txt']);
        
        newfirst = importdata([resultspath name{i} 'Seed' num2str(seed(ii)) 'FirstResponder' num2str(val(v)) '.txt']);
        %secondresponder = importdata([name{i} 'Seed' num2str(seed(ii)) 'SecondResponder.txt']);
        
        memoffirst = ismember(newfirst, origfirst);
        memofsecond = ismember(newfirst, origsecond);
        
        Data.name{count, 1} = name{i};
        Data.seed(count,1) = seed(ii);
        Data.percentfirst(count,1) = sum(memoffirst)/length(origfirst);
        Data.percentsecond(count,1) = sum(memofsecond)/length(origfirst);
        count=count+1;
        
    end
end

filename = [mainpath savename datestr(datetime('today'))];
T = struct2table(Data);
writetable(T,[filename '.xlsx']);
%