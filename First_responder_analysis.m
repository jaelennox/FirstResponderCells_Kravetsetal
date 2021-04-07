
%% THIS PROGRAM CALCULATES TIME OF THE 1-ST PHASE RESPONCE OF THE CALCIUM SIGNAL OF THE ISLET AND EACH INDIVIDUAL CELL
%% INPUT PARAMETERS: NUMBER OF CELLS, TIME INTERVAL, CALCIUM TIMELAPSE DATA
%% CREDIT: VIRA KRAVETS AUG, 2019
close all
clear all
clc

random = 0;

%% 1. INPUT PARAMETERS
seed = {'1'};
name = 'WT100';
checkcells = 0;
val = [100];
files = {'WT_G2to11_Coup120Seed1/'};

for v = 1:length(val)
    for ii = 1:length(files)
        filepath = [files{ii}];
        calciumT = importdata([filepath '/calcium.txt']);    %update this path for your folder
        
        RandomVarsT=importdata([filepath '/RandomVars.txt']);
        
        if checkcells
            ZeroCoupCell = importdata([filepath '/ZeroCoupCell.txt']); %use this line for "PostAblation" analysis
            ZeroCoupCell = ZeroCoupCell+1;
            allcells = ones(1,1000);
            allcells(ZeroCoupCell) = 0;
            nonZero = find(allcells==1);
            if random
                newfile = regexprep(filepath,'_random','');
                ZeroCoupCell = importdata([newfile '/ZeroCoupCell.txt']); %use this line for "PostAblation" analysis
                ZeroCoupCell = ZeroCoupCell+1;
            end
        end
        
        numcells=1000;                                       % enter number of cells
        
        figure
        plot(calciumT)
        st=200;                                              % starting frame (low glucose before "swap")
        ed=720;                                              % ending frame (high glucose)
        
        figure
        plot(calciumT(st:ed,:))                              % plot Ca as is
        [Mx,IndMx]=maxk(calciumT(st:ed,:),1,1);              % returns max values of the Ca intensity for each cell, and corresponding index (time point)
        [Mn,IndMn]=mink(calciumT(st:ed,:),1,1);
        
        % 3. MAKING THE REFERENCE SIGNAL TO COMPARE THE SIGNAL OF INDIVIDUAL CELL'S Ca WITH THIS REFERENCE
        
        for i=st:ed                                          % itterative index for time
            currenttime = calciumT(i,:);                         % picking the row of the calcium table, corresponding to current timepoint, i
            MeanIslet(i-(st-1)) = mean(currenttime);             % reference signal. Index (i-(st-1)) is here to account for times when st is not 0, otherwise indexing is wrong
            
        end
        figure
        plot(MeanIslet)                                      % plotting reference signal (Islet-average)
        
        % 4. NORMALIZING Ca OF ISLET_AVERAGE AND EACH INDIVIDUAL CELL TO BE BETWEEN [0:1]
        % and OBTAINING CROSS-CORRELATION OF THE REFERENCE SIGNAL (MEANISLET)WITH EACH INDIVIDUAL CELL
        
        MeanIsletN=((MeanIslet-MeanIslet(1,1))./(mean(Mx)-MeanIslet(1,1)));  % normalizing islet-average to be between [0:1]
        getdistancefromhalf = MeanIsletN-0.5;
        HHtime = min(find(getdistancefromhalf>0));
        %[HHval,HHtime] = min(abs(MeanIsletN-0.5));           % HHtime - time at which Ca elevation of the Islet-Average reaches it's half-height; HHval - not important (equals to [normalized Ca intensity - 0.5])
        
        for k=1:numcells                                     % itterative index for cells
            currentcell = calciumT(st:ed,k);                 % picking column of the calcium table corresponding to current cell, j
            currentcellN =((currentcell-Mn(k))./(Mx(k)-Mn(k)));   % normalizing each cell to be between [0:1]
            calciumN(1:ed-st+1,k)=currentcellN;              % writing each normalized cell into array
            
            %     [c(:,k) lags]=xcov(currentcellN,MeanIsletN,'coeff');    % cross-covariance  measures the similarity between currentcell and shifted (lagged) copies of MeanIslet as a function of the lag.
            %     [maxCV maxCL]=max(c);                                   % maxCV - maximum value of xcov; maxCL - index (lag) of the maxCV
            getdistancefromhalf = currentcellN-0.5;
            cHHtime(:,k) = min(find(getdistancefromhalf>0));
            %[cHHval,cHHtime(:,k)] = min(abs(currentcellN-0.5));  % cHHtime - time at which Ca elevation of the k-th cell reaches it's half-height; HHval - not important
        end
        
        figure
        plot(calciumN)                                       % plotting normalized [0:1] Ca timelapses
        
        [FirstVal,FirstInd]=mink(cHHtime,1000,2);            % FirstVal - HHElevT of the 1000 (or other number) 1st-responders; FirstInd - cell numbers of the first 1000 (or other number) 1st-responders
        FirstIndCol=FirstInd';                               % converting cell numbers of the 1st responding cells to column vector
        Tresp=cHHtime';                                      % converting time at which Ca elevation of the k-th cell reaches it's half-height to column vector
        
        %get AUC for the increase
        AUC(ii) = trapz(MeanIslet(HHtime:HHtime+300));
        AvgTresp(ii) = HHtime; % Period of the main component
        %
        % 6. PLOTTING gKatp and kGlyc FOR THE 1st RESPONDERS AND FOR THE ISLET-AVERAGE.
        gKatp=RandomVarsT(:,1);              % pulling Katp (potassium channel conductance) from RandomVars's column 1;
        gKatpMean(ii)=mean(gKatp);
        gKatpStd(ii)=std(gKatp);
        
        gCoup=RandomVarsT(:,2);              % pulling gCoup (coupling conductance) from RandomVars's column 2;
        gCoupMean(ii)=mean(gCoup);
        gCoupStd(ii)=std(gCoup);
        
        kGlyc=RandomVarsT(:,10);              % pulling kGlyc (glycolysis rate) from RandomVars's column 10;
        kGlycMean(ii)=mean(kGlyc);
        kGlycStd(ii)=std(kGlyc);
        
        %%%%%Locating first responders and their corresponding metabolic and
        %%%%%electrical parameters:
        jj = 1;
        for j=FirstInd(1:val(v))                 % here FirstInd(1:x) includes x=max number of first-responder cells to average over
            gKatpCellsF(jj,ii)=gKatp(j);         % locating gKatp corresponding to the first x fastest-responding cells
            gCoupCellsF(jj,ii)=gCoup(j);
            kGlycCellsF(jj,ii)=kGlyc(j);
            TrespF(jj,ii)=Tresp(j);
            jj = jj+1;% locating kGlyc corresponding to the first x fastest-responding cells
        end
        
        jj = 1;
        for j=FirstInd(val(v)+1:numcells)                 % here FirstInd(1:x) includes x=max number of first-responder cells to average over
            gKatpCellsNonF(jj,ii)=gKatp(j);         % locating gKatp corresponding to the first x fastest-responding cells
            gCoupCellsNonF(jj,ii)=gCoup(j);
            kGlycCellsNonF(jj,ii)=kGlyc(j);
            TrespNonF(jj,ii)=Tresp(j);
            jj = jj+1;% locating kGlyc corresponding to the first x fastest-responding cells
        end
        
        if checkcells
            jj = 1;
            for j=ZeroCoupCell'% here FirstInd(1:x) includes x=max number of first-responder cells to average over
                gKatpCellsZ(jj,ii)=gKatp(j);         % locating gKatp corresponding to the first x fastest-responding cells
                gCoupCellsZ(jj,ii)=gCoup(j);
                kGlycCellsZ(jj,ii)=kGlyc(j);
                TrespZ(jj,ii)=Tresp(j);
                jj = jj+1;% locating kGlyc corresponding to the first x fastest-responding cells
            end
            
            jj = 1;
            for j=nonZero                % here FirstInd(1:x) includes x=max number of first-responder cells to average over
                gKatpCellsNonZ(jj,ii)=gKatp(j);         % locating gKatp corresponding to the first x fastest-responding cells
                gCoupCellsNonZ(jj,ii)=gCoup(j);
                kGlycCellsNonZ(jj,ii)=kGlyc(j);
                TrespNonZ(jj,ii)=Tresp(j);
                jj = jj+1;% locating kGlyc corresponding to the first x fastest-responding cells
            end
        end
        
        % here FirstInd(1:x) includes x=max number of first-responder cells to average over       % locating kGlyc corresponding to the first x fastest-responding cell
        
        gKatpCellsNZ=nonzeros(gKatpCellsF);    % removing all the cells not belonging to fastest [1:x]
        FirstgKatpMean=mean(gKatpCellsNZ);    % finding avreage gKatp for the first [1:x] cells
        
        gCoupCellsNZ=nonzeros(gCoupCellsF);    % removing all the cells not belonging to fastest [1:x]
        FirstgCoupMean=mean(gCoupCellsNZ);   % finding avreage gCoup for the first [1:x] cells
        
        kGlycCellsNZ=nonzeros(kGlycCellsF);    % removing all the cells not belonging to fastest [1:x]
        FirstgkGlycMean=mean(kGlycCellsNZ);   % finding avreage kGlyc for the first [1:x] cells
        
        first_parameters(ii,:) = [FirstgKatpMean,  FirstgCoupMean, FirstgkGlycMean];
        
        %%%%%Locating last responders and their corresponding metabolic and
        %%%%%electrical parameters:
        jj=1;
        for j=FirstInd(numcells-val(v)+1:numcells)                 % here FirstInd(x:numcells) includes x= number of last-responder cells to average over
            gKatpCellsL(jj,ii)=gKatp(j);         % locating gKatp corresponding to the first x last-responding cells
            gCoupCellsL(jj,ii)=gCoup(j);
            kGlycCellsL(jj,ii)=kGlyc(j);
            TrespL(jj,ii)=Tresp(j);
            jj=jj+1;% locating kGlyc corresponding to the first x last-responding cells
        end
        
        jj=1;
        for j=FirstInd(1:numcells-val(v))                 % here FirstInd(x:numcells) includes x= number of last-responder cells to average over
            gKatpCellsNonL(jj,ii)=gKatp(j);         % locating gKatp corresponding to the first x last-responding cells
            gCoupCellsNonL(jj,ii)=gCoup(j);
            kGlycCellsNonL(jj,ii)=kGlyc(j);
            TrespNonL(jj,ii)=Tresp(j);
            jj=jj+1;% locating kGlyc corresponding to the first x last-responding cells
        end
        
        gKatpCellsNZ=nonzeros(gKatpCellsL);    % removing all the cells not belonging to slowest [1:x]
        LastKatpMean=mean(gKatpCellsNZ);    % finding avreage gKatp for the last [1:x] cells
        
        gCoupCellsNZ=nonzeros(gCoupCellsL);    % removing all the cells not belonging to slowest [1:x]
        LastgCoupMean=mean(gCoupCellsNZ);    % finding avreage gKatp for the last [1:x] cells
        
        kGlycCellsNZ=nonzeros(kGlycCellsL);    % removing all the cells not belonging to slowest [1:x]
        LastgkGlycMean=mean(kGlycCellsNZ);   % finding avreage kGlyc for the last [1:x] cells
        
        figure
        scatter(Tresp,gKatp);
        title('gKatp, Katp conductance')
        figure
        scatter(Tresp,kGlyc);
        title('GLYCv, Glycolysis Rate')
        
        %%%%%%%%%
        FirstRespCa=calciumT(:,FirstInd(1));          %finding the ca timelapse corresponding to 1st-responder
        LastRespCa=calciumT(:,FirstInd(end));         %finding the ca timelapse corresponding to last-responder
        figure
        plot(FirstRespCa);
        title('Calcium for 1st, mean, and last -responders')
        hold on
        plot(mean(calciumT,2));
        hold on
        plot(LastRespCa)
        
    end
    
    %%means of first and last
    gKatpCellsFirst_mean=mean(gKatpCellsF)';         % locating gKatp corresponding to the first x fastest-responding cells
    gCoupCellsFirst_mean=mean(gCoupCellsF)';
    kGlycCellsFirst_mean=mean(kGlycCellsF)';
    TrespCellsFirst_mean=mean(TrespF)';
    gKatpCellsLast_mean=mean(gKatpCellsL)';       % locating gKatp corresponding to the first x fastest-responding cells
    gCoupCellsLast_mean=mean(gCoupCellsL)';
    kGlycCellsLast_mean=mean(kGlycCellsL)';
    TrespCellsLast_mean=mean(TrespL)';
    
    if checkcells
        gKatpCellsZero_mean=mean(gKatpCellsZ)';         % locating gKatp corresponding to the first x fastest-responding cells
        gCoupCellsZero_mean=mean(gCoupCellsZ)';
        kGlycCellsZero_mean=mean(kGlycCellsZ)';
        TrespCellsZero_mean=mean(TrespZ)';
    end
    
    %%means of nonfirst and last
    gKatpCellsFirstNon_mean=mean(gKatpCellsNonF)';         % locating gKatp corresponding to the first x fastest-responding cells
    gCoupCellsFirstNon_mean=mean(gCoupCellsNonF)';
    kGlycCellsFirstNon_mean=mean(kGlycCellsNonF)';
    TrespCellsFirstNon_mean=mean(TrespNonF)';
    gKatpCellsLastNon_mean=mean(gKatpCellsNonL)';       % locating gKatp corresponding to the first x fastest-responding cells
    gCoupCellsLastNon_mean=mean(gCoupCellsNonL)';
    kGlycCellsLastNon_mean=mean(kGlycCellsNonL)';
    TrespCellsLastNon_mean=mean(TrespNonL)';
    
    randlabels = {'gKATP', 'gCoup', 'kGlyc'};
    
    allrandfirst = [gKatpCellsF, gCoupCellsF, kGlycCellsF];
    
    m = 0;
    for n = 1:length(randlabels)
        allrandfirststruct.([randlabels{n} seed{1}]) = allrandfirst(:,length(seed)*m+1);
        m=m+1;
    end
    
    MeansLabels = {'IsletAvggKATP', 'IsletAvgCoup', 'IsletAvgkGlc','IsletStdgKATP', 'IsletStdCoup', 'IsletStdkGlc',...
        'AUCmeanislet', 'IsletAvgTresp', 'gKatpCellsFirst_mean', 'gCoupCellsFirst_mean', 'kGlycCellsFirst_mean', 'TrespCellsFirst_mean',...
        'gKatpCellsLast_mean', 'gCoupCellsLast_mean', 'kGlycCellsLast_mean', 'TrespCellsLast_mean', 'DeltaTresp'};
    Means = [gKatpMean', gCoupMean', kGlycMean', gKatpStd', gCoupStd', kGlycStd', ...
        AUC'/10*1000, AvgTresp'/10, gKatpCellsFirst_mean, gCoupCellsFirst_mean, kGlycCellsFirst_mean, TrespCellsFirst_mean/10,...
        gKatpCellsLast_mean, gCoupCellsLast_mean, kGlycCellsLast_mean, TrespCellsLast_mean/10, (TrespCellsLast_mean-TrespCellsFirst_mean)/10];
    
    MeansLabelsNon = {'gKatpCellsFirstNon_mean', 'gCoupCellsFirstNon_mean', 'kGlycCellsFirstNon_mean','TrespCellsFirstNon_mean',...
        'gKatpCellsLastNon_mean', 'gCoupCellsLastNon_mean', 'kGlycCellsLastNon_mean', 'TrespCellsLastNon_mean'};
    MeansNon = [gKatpCellsFirstNon_mean, gCoupCellsFirstNon_mean, kGlycCellsFirstNon_mean, TrespCellsFirstNon_mean/10,......
        gKatpCellsLastNon_mean, gCoupCellsLastNon_mean, kGlycCellsLastNon_mean, TrespCellsLastNon_mean/10];
    
    MeanLabels1 = [MeansLabels, MeansLabelsNon];
    Means1 = [Means, MeansNon];
    
    for m = 1:length(Means1)
        allmeans.(MeanLabels1{m}) = Means1(:,m);
    end
    
    %     filename = [resultspath name 'summarydata' datestr(datetime('today'))];
    %     T = struct2table(allmeans);
    %     writetable(T,[filename '.xlsx']);
    
    % save([resultspath 'workspace_first_ind_parameters_' name num2str(val(v))])
end