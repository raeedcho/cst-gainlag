function [Data,miss4] = errorCursorSaveFix(Data)
%ERRORCURSORSAVEFIX Recover missing samples from successful CST trials
%   Fix for data collected with MonkeyHost version 09.28.2017.
%   Missing samples were dropped during sending to Host; these samples were
%   still used during online CST control.
%   Up to three consecutive missing samples can be recovered.
%   Note: This function only interpolates for successful trials.
%
%   Arguments
%   ---------
%   Data : 1xN struct array
%       Data structure.
%
%   Returns
%   ---------
%   Data : 1xN struct array
%       Data structure with corrected CST data.
%   miss4 : vector
%       Trials with >=4 samples still missing.
%
%   Basic Usage::
%       >> [Data,miss4] = errorCursorSaveFix(Data)
%
%   2018 Nicole McClain

%% Extract successful CST trials
trialTypes = arrayfun(@(data) {data.Overview.trialName}, Data);
cstTrials = ~cellfun(@isempty,regexpi(trialTypes,'CST'));
replayTrials = ~cellfun(@isempty,regexpi(trialTypes,'CST Trial Replay'));
success = arrayfun(@(data) data.Overview.trialStatus == 1, Data);
if ~sum(cstTrials & success & ~replayTrials)
    error('No successful CST trials in this session.')
end
%Extract successful CST trial data
trls = find(cstTrials & success & ~replayTrials); %find indices
cstData = Data(trls);

%% Process CST trials
miss4 = false(size(trls)); %set up to check for missing samples
for i = 1:length(cstData)
    %% Extract known trial information
    ref = cstData(i).Parameters.ForceParameters.refValue;
    l = cstData(i).Parameters.ForceParameters.initialLambda;
    pspX = cstData(i).TrialData.Marker.rawPositions(:,2); %x axis input
    stateNames = [cstData(i).Parameters.StateTable.stateName];
    cstState = find(strcmp(stateNames,'Control System'));
    stateTrans = cstData(i).TrialData.stateTransitions;
    %(State time)-1 because Samples occur at mod(trialClock-1,10) == 0
    cstStartT = double(stateTrans(2,stateTrans(1,:)==cstState))-1;
    erc = cstData(i).TrialData.Marker.errorCursor;
    %Handle translation issues (e.g. double-sized or multiple z-values)
    erc = erc(erc(:,3)==erc(1,3),:);
    [x,t,u] = deal(erc(:,1),erc(:,4),erc(:,6));
    idx = (t/10-ceil(cstStartT/10))+1; %sample index
    uPSP = [pspX; 0]-ref; %possible input values
    %Check if samples start at xFull(1) == ref
    if xor(idx(1)==1,x(1)==ref)
        error('Indexing mismatch on trial %i',trls(i));
    end
    
    %% Preallocate CST data
    [xFull,tFull,uFull] = deal(NaN(600,1));
    [xFull(idx),tFull(idx),uFull(idx)] = deal(x,t,u);
    
    %% Calculate missing first sample
    if (idx(1)-1) == 1 %one missing sample
        uPoss = [uPSP(1:find(uPSP == u(1))); 0];
        if length(uPoss)==1; uPoss = uPSP; end
        xPoss = 0+(exp(l*0.01)-1)*uPoss;
        if ~sum(xPoss==x(1)) %check for rounding issues
            error('No matching output for trial %i',trls(i));
        end
        xFull(1) = 0;
        uFull(1) = uPoss(find(xPoss==x(1),1));
        tFull(1) = ceil(cstStartT/10)*10;
    elseif (idx(1)-1) == 2 %two missing samples
        uPoss = [uPSP(1:find(uPSP == u(1))); 0];
        if length(uPoss)==1; uPoss = uPSP; end
        xPoss = 0+(exp(l*0.01)-1)*uPoss;
        [xPoss1,uPossG] = meshgrid(xPoss,uPoss);
        xPoss2 = (exp(l*0.01)*xPoss1+(exp(l*0.01)-1)*uPossG);
        if ~sum(xPoss2==x(1)) %check for rounding issues
            error('No matching output for trial %i',trls(i))
        end
        xFull(2) = xPoss1(find(xPoss2==x(1),1));
        uFull(2) = uPossG(find(xPoss2==x(1),1));
        xFull(1) = 0;
        uFull(1) = uPoss(find(xPoss==xFull(2),1));
        tFull(1:2) = ceil(cstStartT/10)*10+([0 10]);
    elseif (idx(1)-1) == 3 %three missing samples
        uPoss = [uPSP(1:find(uPSP == u(1))); 0];
        if length(uPoss)==1; uPoss = uPSP; end
        xPoss = 0+(exp(l*0.01)-1)*uPoss;
        [xPoss1,uPossG] = meshgrid(xPoss,uPoss);
        xPoss2 = (exp(l*0.01)*xPoss1+(exp(l*0.01)-1)*uPossG);
        [xPoss2G,uPossG3] = meshgrid(xPoss2,uPoss);
        xPoss3 = (exp(l*0.01)*xPoss2G+(exp(l*0.01)-1)*uPossG3);
        if ~sum(sum(sum(xPoss3==x(1)))) %check for rounding issues
            error('No matching output for trial %i',trls(i))
        end
        xFull(3) = xPoss2G(find(xPoss3==x(1),1));
        uFull(3) = uPossG3(find(xPoss3==x(1),1));
        xFull(2) = xPoss1(find(xPoss2==xFull(3),1));
        uFull(2) = uPossG(find(xPoss2==xFull(3),1));
        xFull(1) = 0;
        uFull(1) = uPossG(find(xPoss1==xFull(2),1));
        tFull(1:3) = ceil(cstStartT/10)*10+([0 10 20]);
    elseif (idx(1)-1) > 3 %too many missing samples to fix
        miss4(i) = true;
    end
    
    %% Fill in missing samples
    nMiss = diff(idx)-1; %number of missing samples
    for j = 1:length(nMiss)
        %check potential output for all
        if nMiss(j) == 1 %simple case, one missing sample
            uPoss = [uPSP(find(uPSP == u(j)):find(uPSP == u(j+1))); 0];
            if length(uPoss)==1; uPoss = uPSP; end
            xNext = exp(l*0.01)*x(j)+(exp(l*0.01)-1)*u(j);
            xPoss = (exp(l*0.01)*xNext+(exp(l*0.01)-1)*uPoss);
            xFull(idx(j)+1) = xNext;
            if ~sum(xPoss==x(j+1)) %check for rounding issues
                error('No matching output for trial %i',trls(i))
            end
            uFull(idx(j)+1) = uPoss(find(xPoss==x(j+1),1));
            tFull(idx(j)+1) = t(j)+10;
        elseif nMiss(j) == 2 %two missing samples
            uPoss = [uPSP(find(uPSP == u(j)):find(uPSP == u(j+1))); 0];
            if length(uPoss)==1; uPoss = uPSP; end
            xNext = exp(l*0.01)*x(j)+(exp(l*0.01)-1)*u(j);
            xPoss1 = (exp(l*0.01)*xNext+(exp(l*0.01)-1)*uPoss);
            [xPoss1,uPossG] = meshgrid(xPoss1,uPoss);
            xPoss2 = (exp(l*0.01)*xPoss1+(exp(l*0.01)-1)*uPossG);
            xFull(idx(j)+1) = xNext;
            if ~sum(sum(xPoss2==x(j+1))) %check for rounding issues
                error('No matching output for trial %i',trls(i));
            end
            uFull(idx(j)+2) = uPossG(find(xPoss2==x(j+1),1));
            xFull(idx(j)+2) = xPoss1(find(xPoss2==x(j+1),1));
            uFull(idx(j)+1) = uPossG(find(xPoss1==xFull(idx(j)+2),1));
            tFull(idx(j)+(1:2)) = t(j)+10*(1:2);
        elseif nMiss(j) == 3
            uPoss = [uPSP(find(uPSP == u(j)):find(uPSP == u(j+1))); 0];
            if length(uPoss)==1; uPoss = uPSP; end
            xNext = exp(l*0.01)*x(j)+(exp(l*0.01)-1)*u(j);
            xPoss1 = (exp(l*0.01)*xNext+(exp(l*0.01)-1)*uPoss);
            [xPoss1,uPossG] = meshgrid(xPoss1,uPoss);
            xPoss2 = (exp(l*0.01)*xPoss1+(exp(l*0.01)-1)*uPossG);
            [xPoss2G,uPossG3] = meshgrid(xPoss2,uPoss);
            xPoss3 = (exp(l*0.01)*xPoss2G+(exp(l*0.01)-1)*uPossG3);
            xFull(idx(j)+1) = xNext;
            if ~sum(sum(sum(xPoss3==x(j+1)))) %check for rounding issues
                error('No matching output for trial %i',trls(i))
            end
            xFull(idx(j)+3) = xPoss2G(find(xPoss3==x(j+1),1));
            uFull(idx(j)+3) = uPossG3(find(xPoss3==x(j+1),1));
            xFull(idx(j)+2) = xPoss1(find(xPoss2==xFull(idx(j)+3),1));
            uFull(idx(j)+2) = uPossG(find(xPoss2==xFull(idx(j)+3),1));
            uFull(idx(j)+1) = uPossG(find(xPoss1==xFull(idx(j)+2),1));
            tFull(idx(j)+(1:3)) = t(j)+10*(1:3);
        elseif nMiss(j) > 3 %too many missing samples to fix
            miss4(i) = true;
        end
    end
    
    %% Check for additional missing samples
    if sum(isnan(xFull)) > 0 && ~miss4(i)
        error('Samples still missing on trial %i',trls(i))
    end
    
    %% Create full error cursor data array
    newErc = [xFull,ones(600,1)*erc(1,2:3),tFull,ones(600,1)*l,uFull];
    cstData(i).TrialData.Marker.errorCursor = newErc;
end

%% Report trials with missing samples
miss4 = trls(miss4);
if ~isempty(miss4)
    warning('More than three samples missing on trials\n%s',...
        sprintf('%i ',miss4))
end

%% Add processed data back into structure
Data(trls) = cstData;

end