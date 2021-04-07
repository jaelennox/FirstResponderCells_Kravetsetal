
%% THIS PROGRAM CALCULATES TIME OF THE 1-ST PHASE RESPONCE OF THE CALCIUM SIGNAL OF THE ISLET AND EACH INDIVIDUAL CELL
%% INPUT PARAMETERS: NUMBER OF CELLS, TIME INTERVAL, CALCIUM TIMELAPSE DATA AS TABLE
%% CREDIT: VIRA KRAVETS AUG, 2019
close all
clear all
clc


addpath('C:\Users\dwuletj\Documents\GitHub\UniversalCode\')
addpath('/Users/jdwulet/Documents/GitHub/UniversalCode/')

%% 1. INPUT PARAMETERS
mainpath = '/Volumes/Seagate Backup Plus Drive/ForViraNewResponders/';
resultspath = '/Users/jdwulet/Google Drive/LabStuff/Projects/FirstResponder/';

% name = 'uncoupled_perc10_v2';
% files = {'Uncouple_FirstResp_Seed1_perc10', 'Uncouple_FirstResp_Seed2_perc10', 'Uncouple_FirstResp_Seed3_perc10', 'Uncouple_FirstResp_Seed4_perc10', 'Uncouple_FirstResp_Seed5_perc10'};
% 
% name = 'uncoupled_perc05_v2';
% files = {'Uncouple_FirstResp_Seed1_perc05', 'Uncouple_FirstResp_Seed2_perc05', 'Uncouple_FirstResp_Seed3_perc05', 'Uncouple_FirstResp_Seed4_perc05', 'Uncouple_FirstResp_Seed5_perc05'};
% name = 'uncoupled_perc1_v2';
% files = {'Uncouple_FirstResp_Seed1_perc1', 'Uncouple_FirstResp_Seed2_perc1', 'Uncouple_FirstResp_Seed3_perc1', 'Uncouple_FirstResp_Seed4_perc1', 'Uncouple_FirstResp_Seed5_perc1'};
% name = 'uncoupled_perc20_v2';
% files = {'Uncouple_FirstResp_Seed1_perc20', 'Uncouple_FirstResp_Seed2_perc20', 'Uncouple_FirstResp_Seed3_perc20', 'Uncouple_FirstResp_Seed4_perc20', 'Uncouple_FirstResp_Seed5_perc20'};
% name = 'uncoupled_perc30_v2';
% files = {'Uncouple_FirstResp_Seed1_perc30', 'Uncouple_FirstResp_Seed2_perc30', 'Uncouple_FirstResp_Seed3_perc30', 'Uncouple_FirstResp_Seed4_perc30',  'Uncouple_FirstResp_Seed5_perc30'};
% name = 'uncoupled_perc50_v2';
% files = {'Uncouple_FirstResp_Seed1_perc50', 'Uncouple_FirstResp_Seed2_perc50', 'Uncouple_FirstResp_Seed3_perc50', 'Uncouple_FirstResp_Seed4_perc50', 'Uncouple_FirstResp_Seed5_perc50'};
% 
% name = 'uncoupled_random_perc05_v2';
% files = {'AblationSeed1_perc05', 'AblationSeed2_perc05', 'AblationSeed3_perc05', 'AblationSeed4_perc05', 'AblationSeed5_perc05'};
% name = 'uncoupled_random_perc1_v2';
% files = {'AblationSeed1_perc1', 'AblationSeed2_perc1', 'AblationSeed3_perc1', 'AblationSeed4_perc1',  'AblationSeed5_perc1'};
% name = 'uncoupled_random_perc10_v2';
% files = {'AblationSeed1_perc10', 'AblationSeed2_perc10', 'AblationSeed3_perc10', 'AblationSeed4_perc10', 'AblationSeed5_perc10'};
% name = 'uncoupled_random_perc20_v2';
% files = {'AblationSeed1_perc20', 'AblationSeed2_perc20', 'AblationSeed3_perc20', 'AblationSeed4_perc20', 'AblationSeed5_perc20'};
% name = 'uncoupled_random_perc30_v2';
% files = {'AblationSeed1_perc30', 'AblationSeed2_perc30', 'AblationSeed3_perc30', 'AblationSeed4_perc30', 'AblationSeed5_perc30'};
% name = 'uncoupled_random_perc50_v2';
% files = {'AblationSeed1_perc50', 'AblationSeed2_perc50', 'AblationSeed3_perc50', 'AblationSeed4_perc50', 'AblationSeed5_perc50'};

%%change XYZ
name = 'uncoupled_perc10_changeXYZ';
files = {'Uncouple_FirstResp_Seed1_perc10_changeXYZ', 'Uncouple_FirstResp_Seed2_perc10_changeXYZ', 'Uncouple_FirstResp_Seed3_perc10_changeXYZ', 'Uncouple_FirstResp_Seed4_perc10_changeXYZ', 'Uncouple_FirstResp_Seed5_perc10_changeXYZ'};

name = 'uncoupled_perc10_changeXYZ_newsphere';
files = {'Uncouple_FirstResp_Seed1_perc10_changeXYZ_newsphere', 'Uncouple_FirstResp_Seed2_perc10_changeXYZ_newsphere', 'Uncouple_FirstResp_Seed3_perc10_changeXYZ_newsphere', 'Uncouple_FirstResp_Seed4_perc10_changeXYZ_newsphere', 'Uncouple_FirstResp_Seed5_perc10_changeXYZ_newsphere'}


seed = {'1', '2', '3', '4', '5'};
numfirstresponder = 100;

%% 1. INPUT PARAMETERS
for ii = 1:length(files)
    
    filepath = [mainpath files{ii}];
    %filepath = '/Users/vira/Desktop/Simulations_WorkFromHome11_20_19/Islet1_leadingVKcorrectSeed1';
    
    calciumT = importdata([filepath '/calcium.txt']);
    %calciumT = table2array(calcium);
    RandomVarsT=importdata([filepath '/RandomVars.txt']);
    %RandomVarsT = table2array(RandomVars);
    ZeroCoupCell = importdata([filepath '/ZeroCoupCell.txt']);
    ZeroCoupCell = ZeroCoupCell+1;                   % takes into account difference in Matlab and C++;
    
    ZeroCoupCellSorted = sort(ZeroCoupCell);         % sorts ZeroCoupCell (cell numbers) in ascending order
    
    numcells=size(calciumT,2);                       % number of cells determined automatially
    
    figure(1)
    plot(calciumT)
    title('Calcium all cells')
    
    st=200;                                          % starting frame (2->11 mM swap time)
    ed=720;                                          % ending frame (Ca platoe time, or if platoe is different for different cells, ed = max time frame)
    
    stToed = ed-st;
    timesteps=size(calciumT,1);
    numuncoupled = length(ZeroCoupCell);
    UncoupledCa=zeros(timesteps, numuncoupled);
    
    count = 1;
    for i=1:numuncoupled
        UncoupledCa(:,count)=calciumT(:,ZeroCoupCell(i));  % timelapse for uncoupled cell
        count=count+1;
    end
    
    figure(2)
    plot(UncoupledCa);
    title('Calcium Uncoupled cells')
    
    % obtaining cell numbers of the coupled cells
    CoupCells = [];
    numcoupled = numcells-numuncoupled;   % for ex: #coupled = 1000-101 = 899
    for i=1:numcells                      % for ex: for i=1:1000
        if ~ismember(i,ZeroCoupCell)      % if i is NOT (~ means NOT) a member of the ZeroCoupCell array
            CoupCells = [CoupCells; i];   % write the cell number into the CoupCells array
        end
    end
    
    % obtaining corresponding Ca for the coupled cells
    CoupledCa = [];
    count = 1;
    for i=1:numcoupled
        CoupledCa(:,count)=calciumT(:,CoupCells(i)); % Calcium timecourse for the coupled cells
        count=count+1;
    end
    
    figure(3)
    plot(CoupledCa);
    title('Calcium Coupled cells, whole T range')
    
    
    % 3. MAKING THE REFERENCE SIGNAL TO COMPARE THE SIGNAL OF INDIVIDUAL CELL'S Ca WITH THIS REFERENCE
    
    for i=st:ed                                       % itterative index for time
        currenttime = CoupledCa(i,:);                     % picking the row of the calcium table, corresponding to current timepoint, i
        MeanIslet(i-(st-1)) = mean(currenttime);          % reference signal (based only on coupled cells). Index (i-(st-1)) is here to account for times when st is not 0, otherwise indexing is wrong
    end
    
    figure(5)
    plot(MeanIslet)
    title('Calcium Mean [Coupled cells]')
    
    % 4. NORMALIZING Ca OF ISLET_AVERAGE AND EACH INDIVIDUAL CELL TO BE BETWEEN [0:1]
    
    % finding max and min of each cell's Ca for normalization
    [Mx,IndMx]=maxk(CoupledCa(st:ed,:),1,1);             % returns 1 max values of the 1st dimension (Ca intensity) for st:ed times for all cells, and corresponding index (time point)
    [Mn,IndMn]=mink(CoupledCa(st:ed,:),1,1);
    
    MeanIsletN=((MeanIslet-MeanIslet(1,1))./(mean(Mx)-MeanIslet(1,1)));  % normalizing islet-average to be between [0:1]
    getdistancefromhalf = MeanIsletN-0.5;
    HHtime = min(find(getdistancefromhalf>0));
    %[HHval,HHtime] = min(abs(MeanIsletN(1:400)-0.5));                   % HHtime - time at which Ca elevation of the Islet-Average reaches it's half-height; HHval - normalized Ca intensity - 0.5 = approaches 0 at HHtime.
    % in the line above "abs" insures that half of the Ca curve below 0 is compared to the half of the Ca curve above 0 appropriately.
    
    %%%%%
    %normalizing each cell's Ca
    
    for k=1:numcoupled
        currentcell = CoupledCa(st:ed,k);                        % picking column of the calcium table corresponding to current cell, j
        currentcellN =((currentcell-Mn(k))./(Mx(k)-Mn(k)));      % normalizing each cell to be between [0:1]
        calciumN(1:ed-st+1,k)=currentcellN;                      % writing each normalized cell into array
        getdistancefromhalf = currentcellN-0.5;
        cHHtime(:,k) = min(find(getdistancefromhalf>0));
        %[cHHval,cHHtime(:,k)] = min(abs(currentcellN(1:400)-0.5));     % cHHtime - time at which Ca elevation of the k-th cell reaches it's half-height; HHval - normalized Ca intensity - 0.5 = approaches 0 at cHHtime
        
    end
    
    figure(6)
    plot(calciumN)
    title('Calcium Coupled cells Normalized')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% code under construction %%%%%%%%%%%%%%%%
    % 5. FINDING PLATOE TIME: condition for platoe time is: at platoe time the
    % derivative of Ca becomes negative & Ca intensity drops more than 50%.
    
    %5A Derivative of coupled Ca for finding the time range for normalization
    
    CoupledCadT = [];       % creates empty array
    count = 1;              % counts timepoints
    for i=1:numcoupled
        CoupledCadT(:,count)=diff(calciumN(:,i)); % derivative of Calcium timecourse for the coupled cells
        count=count+1;
    end
    
    figure (4)
    plot(CoupledCadT)
    title('derivative of Calcium Coupled cells, whole T range')
    % 5A. finding times for which derinative is negative and setting them to 1, and
    % if derivative is positive - setting it to 0.
    
    for i=1:numcoupled                                   % itterative index for coupled cells
        currentcelldT = CoupledCadT(1:stToed,i);         % derivative of the Ca for current cell (i)
        [timedTneg,col]=find(currentcelldT<0);           % locating times at which derivative is negative
        for j=1:stToed                                    % j is itterative index for timesteps
            if ismember(j,timedTneg)                     % if j is a memeber of array containing times correcponding to negative derivative
                currentcelldT(j)=1;                          % then write 1 in the currentcelldT for that time
            else currentcelldT(j)=0;                     % otherwise - write 0
            end
        end
        CoupCadTbin(1:stToed,i) = currentcelldT;           % create an array containing all coupled cells, where 1 = neg, 0 = positive Ca derivatives
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% end of code under construction %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [FirstVal,FirstInd]=mink(cHHtime,1000,2);                    % FirstVal - Half Height Elev Time of the 1000 (or other number) 1st-responders; FirstInd - cell numbers of the first 100 1st-responders
    FirstIndCol=FirstInd';
    Tresp=cHHtime';
    
    %get AUC for the increase
    AUC(ii) = trapz(MeanIslet(HHtime:HHtime+300));
    AvgTresp(ii) = HHtime;
    
    
    % 6. PLOTTING gKatp and kGlyc FOR THE 1st RESPONDERS AND FOR THE ISLET-AVERAGE.
    gKatp=RandomVarsT(:,1);
    gKatpCoupled = [];
    count = 1;
    for i=1:numcoupled
        gKatpCoupled(count)=gKatp(CoupCells(i));
        count=count+1;
    end
    gKatpMean=mean(gKatpCoupled);
    
    gCoup=RandomVarsT(:,2); % pulling gCoup (coupling conductance) from RandomVars
    gCoupCoupled = [];
    count = 1;
    for i=1:numcoupled
        gCoupCoupled(count)=gCoup(CoupCells(i));
        count=count+1;
    end
    gCoupMean=mean(gCoupCoupled);
    
    kGlyc=RandomVarsT(:,10);
    kGlycCoupled = [];
    count = 1;
    for i=1:numcoupled
        kGlycCoupled(count)=kGlyc(CoupCells(i));
        count=count+1;
    end
    kGlycMean=mean(kGlycCoupled);
    
    %%%%%Locating first responders and their corresponding metabolic and electrical parameters:
    
    jj = 1;
    for j=FirstInd(1:numfirstresponder)                 % here FirstInd(1:x) includes x=max number of first-responder cells to average over
        gKatpCellsF(jj,ii)=gKatpCoupled(j);         % locating gKatp corresponding to the first x fastest-responding cells
        gCoupCellsF(jj,ii)=gCoupCoupled(j);
        kGlycCellsF(jj,ii)=kGlycCoupled(j);
        TrespF(jj,ii)=Tresp(j);
        jj = jj+1;% locating kGlyc corresponding to the first x fastest-responding cells
    end
    
    %gKatpCellsNZ=nonzeros(gKatpCellsF);    % removing all the cells not belonging to fastest [1:x]
    FirstgKatpMean=mean(gKatpCellsF);    % finding avreage gKatp for the first [1:x] cells
    
    %gCoupCellsNZ=nonzeros(gCoupCellsF);    % removing all the cells not belonging to fastest [1:x]
    FirstgCoupMean=mean(gCoupCellsF);   % finding avreage gCoup for the first [1:x] cells
    
    %kGlycCellsNZ=nonzeros(kGlycCellsF);    % removing all the cells not belonging to fastest [1:x]
    FirstgkGlycMean=mean(kGlycCellsF);   % finding avreage kGlyc for the first [1:x] cells
    
    %%%%%Locating last responders and their corresponding metabolic and
    %%%%%electrical parameters:
    jj=1;
    for j=FirstInd(numcoupled-numfirstresponder+1:numcoupled)                 % here FirstInd(1:x) includes x= number of last-responder cells to average over
        gKatpCellsL(jj,ii)=gKatpCoupled(j);         % locating gKatp corresponding to the first x fastest-responding cells
        gCoupCellsL(jj,ii)=gCoupCoupled(j);
        kGlycCellsL(jj,ii)=kGlycCoupled(j);
        TrespL(jj,ii)=Tresp(j);
        jj=jj+1;% locating kGlyc corresponding to the first x last-responding cells
    end
    %gKatpCellsNZ=nonzeros(gKatpCellsL);    % removing all the cells not belonging to slowest [1:x]
    LastKatpMean=mean(gKatpCellsL);    % finding avreage gKatp for the last [1:x] cells
    
    %gCoupCellsNZ=nonzeros(gCoupCellsL);    % removing all the cells not belonging to slowest [1:x]
    LastgCoupMean=mean(gCoupCellsL);    % finding avreage gKatp for the last [1:x] cells
    
    %kGlycCellsNZ=nonzeros(kGlycCellsL);    % removing all the cells not belonging to slowest [1:x]
    LastgkGlycMean=mean(kGlycCellsL);   % finding avreage kGlyc for the last [1:x] cells
    
    figure(7)
    scatter(Tresp,gKatpCoupled);
    title('gKatp, Katp conductance')
    figure(8)
    scatter(Tresp,kGlycCoupled);
    title('GLYCv, Glycolysis Rate')
    
    
    Firstvalues = CoupCells(FirstInd(1:numfirstresponder))-1;
%     dlmwrite([resultspath name 'Seed' seed{ii} 'FirstResponder' num2str(numuncoupled) '.txt'], Firstvalues);
    try
        Secondvalues = CoupCells(FirstInd(numfirstresponder+1:numfirstresponder*2))-1;
%         dlmwrite([resultspath name 'Seed' seed{ii} 'SecondResponder' num2str(numuncoupled) '.txt'], Secondvalues);
    catch
    end
   % saveAllFigsToPPT([resultspath name 'Seed' seed{ii} datestr(datetime('today'))]);
end

TrespCellsFirst_mean=mean(TrespF)';
TrespCellsLast_mean=mean(TrespL)';
MeansLabels = {'TrespCellsFirst_mean', 'TrespCellsLast_mean','DeltaTresp', 'AUCmeanislet', 'IsletAvgTresp'};
Means = [TrespCellsFirst_mean/10, TrespCellsLast_mean/10, (TrespCellsLast_mean-TrespCellsFirst_mean)/10, AUC'/10*1000, AvgTresp'/10];


for m = 1:length(Means)
    allmeans.(MeansLabels{m}) = Means(:,m);
end

%whichcells_first;
% filename = [resultspath name 'summarydata' datestr(datetime('today'))];
% T = struct2table(allmeans);
% writetable(T,[filename '.xlsx']);
% 
% save([resultspath name 'workspace_first_ind_parameters_' num2str(numuncoupled)])
