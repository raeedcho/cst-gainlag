function td = convertSMILEtoTD(smile_data)
% CONVERTSMILETOTD converts SMILE trial data type to TrialData

    % Loop over smile_data, which is already split up by trial
    td_cell = cell(1,length(smile_data));
    for trialnum = 1:length(smile_data)
        % initialize with all the stuff we want
        td_cell{trialnum} = struct(...
            'monkey',NaN,...
            'task',NaN,...
            'date_time',NaN,...
            'trial_id',NaN,...
            'result',NaN,...
            'bin_size',NaN,...
            'lambda',NaN,...
            'timevec',NaN,...
            'idx_startTime',NaN,...
            'idx_ctHoldTime',NaN,...
            'idx_goCueTime',NaN,...
            'idx_otHoldTime',NaN,...
            'idx_rewardTime',NaN,...
            'idx_failTime',NaN,...
            'idx_endTime',NaN...
            );
    
        td_cell{trialnum} = parse_smile_meta(td_cell{trialnum},smile_data(trialnum));
        td_cell{trialnum} = parse_smile_events(td_cell{trialnum},smile_data(trialnum));
        td_cell{trialnum} = parse_smile_behavior(td_cell{trialnum},smile_data(trialnum));
        td_cell{trialnum} = parse_smile_spikes(td_cell{trialnum},smile_data(trialnum));
    end
    td=cat(2,td_cell{:});
end

function trial = parse_smile_meta(in_trial,smile_trial)
    % Get meta data out of smile trial data structure

    % copy over input trial
    trial = in_trial;
    trial.monkey = smile_trial.Overview.subjectName;
    trial.date_time = datestr(datenum(num2str(smile_trial.Overview.date),'yyyymmddHHMMSS'),'yyyy/mm/dd HH:MM:SS');
    trial.trial_id = sscanf(smile_trial.Overview.trialNumber,'Trial%d');
    trial.task = smile_trial.Overview.trialName;

    % trial result (0-failure, 1-success, 2-error/abort)
    switch smile_trial.Overview.trialStatus
        case 0
            trial.result = 'F';
        case 1
            trial.result = 'R';
        case 2
            trial.result = 'A';
    end

    % get lambda for CST trials
    if contains(smile_trial.Overview.trialName,'CST')
        trial.lambda = smile_trial.Parameters.ForceParameters.initialLambda;
    end
end

function trial = parse_smile_events(in_trial,smile_trial)
    %% parse out events from smile trial data structure
    % assume that we only care about stuff from time zero (which should end up being idx_startTime) to end of trial

    % copy over input trial
    trial = in_trial;

    % figure out timevector for this trial
    state_transitions = smile_trial.TrialData.stateTransitions;
    endTime = state_transitions(2,state_transitions(1,:)==-1);
    trial.timevec = 0:ceil(double(endTime));

    % read out end time index
    trial.idx_endTime = interp1(trial.timevec,1:length(trial.timevec),double(endTime),'nearest');

    % reward
    if smile_trial.Overview.trialStatus==1
        trial.idx_rewardTime = timevec_event_index(smile_trial,'Success',trial.timevec);
    end

    if contains(smile_trial.Overview.trialName,'CST')
        % startTime is probably state 1, but it generally seems to be called 'Reach to Center' or 'Move to Center'
        trial.idx_startTime = timevec_event_index(smile_trial,'Move to Center',trial.timevec);

        % center hold
        trial.idx_ctHoldTime = timevec_event_index(smile_trial,'Hold Center',trial.timevec);

        % go cue
        trial.idx_goCueTime = timevec_event_index(smile_trial,'Control System',trial.timevec);

        % fail
        if smile_trial.Overview.trialStatus==0
            trial.idx_failTime = timevec_event_index(smile_trial,'Failure Display',trial.timevec);
        end
    elseif contains(smile_trial.Overview.trialName,'Center Out')
        % startTime is probably state 1, but it generally seems to be called 'Reach to Center' or 'Move to Center'
        trial.idx_startTime = timevec_event_index(smile_trial,'Reach to Center',trial.timevec);

        % center hold
        trial.idx_ctHoldTime = timevec_event_index(smile_trial,'Center Hold',trial.timevec);

        % go cue
        trial.idx_goCueTime = timevec_event_index(smile_trial,'Reach to Target',trial.timevec);

        % target hold
        trial.idx_otHoldTime = timevec_event_index(smile_trial,'Target Hold',trial.timevec);

        % fail
        if smile_trial.Overview.trialStatus==0
            trial.idx_failTime = timevec_event_index(smile_trial,'Failure (Center)',trial.timevec);
            if isnan(trial.idx_failTime)
                trial.idx_failTime = timevec_event_index(smile_trial,'Target Failure',trial.timevec);
                if isnan(trial.idx_failTime)
                    trial.idx_failTime = timevec_event_index(smile_trial,'Non-Attempt',trial.timevec);
                end
            end
        end
    end
end

function idx = timevec_event_index(smile_trial,event_name,timevec)
    % helper function to figure out index for nearest time in timevec
    state_transitions = smile_trial.TrialData.stateTransitions;
    statenum = find(strcmpi(cat(1,smile_trial.Parameters.StateTable.stateName),event_name));
    query_time = state_transitions(2,state_transitions(1,:)==statenum);
    idx = interp1(timevec,1:length(timevec),double(query_time),'nearest');
    
    if isempty(idx)
        idx = NaN;
    end
end

function trial = parse_smile_behavior(in_trial,smile_trial)
    % parse out behavior (hand position and cursor position) from smile data struct

    % copy over input trial
    trial = in_trial;

    % first get the phasespace marker data
    % phasespace_timevec = smile_trial.TrialData.Marker.rawPositions(:,7);
end

function trial = parse_smile_spikes(in_trial,smile_trial)
    % parse out neural activity from smile data struct
    
    % copy over input trial
    trial = in_trial;
end
