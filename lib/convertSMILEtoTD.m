function trial_data = convertSMILEtoTD(smile_data,params)
% CONVERTSMILETOTD converts SMILE trial data type to TrialData
% Outputs:
%   trial_data - TrialData structure
% Inputs:
%   smile_data - data structure to convert from SMILE type to TrialData
%   params - Parameters struct
%       .bin_size - bin size in seconds to convert to (default: 0.001)
%       .array_alias - cell array for renaming array names.
%           First col is name to replace, second col is replacement name. Rows are different aliases
%       .valid_sort_idx - list of array channel unit indices to keep. Defaults to 1:10
%       .meta - struct containing meta information to add to file (default: emtpy struct)
%           Note that monkey, date, task are extracted from smile_data already.
%
% Written by Raeed Chowdhury Feb 2020

    % check if params was passed in
    if nargin<2
        params = struct();
    end

    % first get unit guide over all of smile_data
    unit_guide = get_smile_unit_guide(smile_data,params);
    
    % fill in missing CST cursor samples
    params.fill_err = false;
%     try
        [smile_data,miss4] = errorCursorSaveFix(smile_data);
%     catch ME
%         params.fill_err = true;
%     end

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
            'ct_location',NaN,...
            'ot_location',NaN,...
            'tgtDir',NaN,...
            'tgtMag',NaN,...
            'idx_startTime',NaN,...
            'idx_ctHoldTime',NaN,...
            'idx_goCueTime',NaN,...
            'idx_otHoldTime',NaN,...
            'idx_cstStartTime',NaN,...
            'idx_cstEndTime',NaN,...
            'idx_rewardTime',NaN,...
            'idx_failTime',NaN,...
            'idx_endTime',NaN...
            );
    
        td_cell{trialnum} = parse_smile_meta(td_cell{trialnum},smile_data(trialnum),params);
        td_cell{trialnum} = parse_smile_events(td_cell{trialnum},smile_data(trialnum),params);
        td_cell{trialnum} = parse_smile_behavior(td_cell{trialnum},smile_data(trialnum),params);
        td_cell{trialnum} = parse_smile_eye_data(td_cell{trialnum},smile_data(trialnum),params);
        td_cell{trialnum} = parse_smile_spikes(td_cell{trialnum},smile_data(trialnum),unit_guide,params);

        td_cell{trialnum} = rmfield(td_cell{trialnum},'timevec');
    end
    % check fields for consistency
    num_fields = cellfun(@(x) length(fieldnames(x)),td_cell);
    if length(unique(num_fields))>1
        warning('something has gone wrong')
    end
    trial_data=cat(2,td_cell{:});
    
    trial_data = reorderTDfields(trial_data);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sub-functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial = parse_smile_meta(in_trial,smile_trial,params)
    % Get meta data out of smile trial data structure

    % default params
    meta = [];
    assignParams(who,params)

    % copy over input trial
    trial = in_trial;
    trial.monkey = smile_trial.Overview.subjectName;
    trial.date_time = datestr(datenum(num2str(smile_trial.Overview.date),'yyyymmddHHMMSS'),'yyyy/mm/dd HH:MM:SS');
    trial.trial_id = sscanf(smile_trial.Overview.trialNumber,'Trial%d');

    % simplify task name
    if startsWith(smile_trial.Overview.trialName,'Center Out')
        trial.task = 'CO';
    elseif startsWith(smile_trial.Overview.trialName,'CST')
        trial.task = 'CST';
    else
        trial.task = smile_trial.Overview.trialName;
    end

    % get center target location
    if startsWith(smile_trial.Overview.trialName,'Center Out')
        ct_reach_state = smile_trial.Parameters.StateTable(strcmpi(vertcat(smile_trial.Parameters.StateTable.stateName),'Reach to Center'));
        center_target_idx = strcmpi(ct_reach_state.StateTargets.names,'co start target');
        trial.ct_location = ct_reach_state.StateTargets.location(center_target_idx,:);
        trial.ct_location(2) = -trial.ct_location(2); % phasespace inverts y-axis for some reason
    elseif startsWith(smile_trial.Overview.trialName,'CST')
        ct_reach_state = smile_trial.Parameters.StateTable(strcmpi(vertcat(smile_trial.Parameters.StateTable.stateName),'Move to Center'));
        center_target_idx = strcmpi(ct_reach_state.StateTargets.names,'start');
        trial.ct_location = ct_reach_state.StateTargets.location(center_target_idx,:);
        trial.ct_location(2) = -trial.ct_location(2); % phasespace inverts y-axis for some reason
    end

    % get reach target location
    if startsWith(smile_trial.Overview.trialName,'Center Out')
        ot_reach_state = smile_trial.Parameters.StateTable(strcmpi(vertcat(smile_trial.Parameters.StateTable.stateName),'Reach to Target'));
        reach_target_idx = strcmpi(ot_reach_state.StateTargets.names,'co reach target');
        trial.ot_location = ot_reach_state.StateTargets.location(reach_target_idx,:);
        trial.ot_location(2) = -trial.ot_location(2); % phasespace inverts y-axis for some reason
        
        trial.tgtDir = atan2d(...
            trial.ot_location(2)-trial.ct_location(2),...
            trial.ot_location(1)-trial.ct_location(1));

        trial.tgtMag = sqrt(sum((trial.ot_location-trial.ct_location).^2));
    end

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

    if ~isempty(meta)
        fnames = fieldnames(meta);
        for fieldnum = 1:length(fnames)
            trial.(fnames{fieldnum}) = meta.(fnames{fieldnum});
        end
    end
end

function trial = parse_smile_events(in_trial,smile_trial,params)
    %% parse out events from smile trial data structure (also add timevec to trial based on end time and bin_size)
    % assume that we only care about stuff from time zero (which should end up being idx_startTime) to end of trial

    % default params
    bin_size = 0.001;
    assignParams(who,params)

    % copy over input trial
    trial = in_trial;

    if ~isfield(trial,'bin_size') || isnan(trial.bin_size)
        trial.bin_size = bin_size;
    end

    % figure out timevector for this trial (in seconds)
    % smile_data specifies everything as milliseconds
    if ~isfield(trial,'timevec') || isnan(trial.timevec)
        state_transitions = smile_trial.TrialData.stateTransitions;
        endTime = state_transitions(2,state_transitions(1,:)==-1);
        end_time_in_s = double(endTime)/1000;
        trial.timevec = 0:bin_size:end_time_in_s;
    end

    % assume last index of the timevec is the end time (might be off by one if SMILE data changes to allow non-int transition times)
    trial.idx_endTime = length(trial.timevec);

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
    % helper function to figure out index for nearest time in timevec (assuming timevec is in seconds)
    state_transitions = smile_trial.TrialData.stateTransitions;
    statenum = find(strcmpi(cat(1,smile_trial.Parameters.StateTable.stateName),event_name));
    query_time = state_transitions(2,state_transitions(1,:)==statenum);
    idx = interp1(timevec,1:length(timevec),double(query_time)/1000,'nearest');
    
    if isempty(idx)
        idx = NaN;
    end
end

function trial = parse_smile_behavior(in_trial,smile_trial,params)
    % parse out behavior (hand position and cursor position) from smile data struct
    % assumes in_trial has timevec based on the endTime found in parse_smile_events

    assert(isfield(in_trial,'timevec'),'Missing temporary timevector in trial data for some reason...')

    % assume we're going for millisecond binning
    bin_size = 0.001;
    fill_err = false;
    assignParams(who,params)
    
    % copy over input trial
    trial = in_trial;

    % first get the phasespace marker data
    phasespace_freq = double(smile_trial.TrialData.Marker.frequency);
    phasespace_data = smile_trial.TrialData.Marker.rawPositions;
    phasespace_mask = [1 -1 1]; %phasespace inverts y axis for some reason
    
    % filtering params
    filter_n = 50;
    kaiser_beta = 20;

    if ~isempty(phasespace_data)
        visual_timevec = phasespace_data(:,7)/1000;
        hand_timevec = phasespace_data(:,6)/1000;
        marker_pos = phasespace_data(:,2:4) .* repmat(phasespace_mask,size(phasespace_data,1),1);

        % resample hand pos data to new timevector
        [temp_resampled,t_resamp] = resample_signals(marker_pos,hand_timevec,struct( ...
            'bin_size',bin_size, ...
            'samprate',phasespace_freq,...
            'filter_n',filter_n,...
            'kaiser_beta',kaiser_beta...
            ));
        trial.hand_pos = interp1(t_resamp,temp_resampled,trial.timevec);
        trial.rel_hand_pos = trial.hand_pos-trial.ct_location;

        % resample cursor pos data to new timevector
        [temp_resampled,t_resamp] = resample_signals(marker_pos(:,1:2),visual_timevec,struct( ...
            'bin_size',bin_size, ...
            'samprate',phasespace_freq,...
            'filter_n',filter_n,...
            'kaiser_beta',kaiser_beta...
            ));
        trial.cursor_pos = interp1(t_resamp,temp_resampled,trial.timevec);
        trial.rel_cursor_pos = trial.cursor_pos-trial.ct_location(1:2);

        % if this is a CST trial, we need to add cursor info during the Control System state (but only if there wasn't an error filling missing samples)
        trial.cst_cursor_command = nan(size(trial.cursor_pos));
        if ~isempty(smile_trial.TrialData.Marker.errorCursor)
            cursor_data = smile_trial.TrialData.Marker.errorCursor;
            cursor_timevec = cursor_data(:,4)/1000; % this is not necessarily regularly spaced but should be after the fix
            cursor_pos = cursor_data(:,1:2);
            cursor_pos(:,2) = -cursor_pos(:,2); %phasespace inverts y-axis for some reason
            
            if fill_err
                warning('Unable to fill missing CST samples for date %s, trying to interpolate samples',trial.date_time)

                % Note: following section is unnecessary after cursor position samples are filled in with Nicole's fix (errorCursorSaveFix)
                % interpolate to uniform sampling rate (assuming an on-average consistent sampling rate)
                dt = (cursor_timevec(end)-cursor_timevec(1))/length(cursor_timevec);
                tGrid = (cursor_timevec(1):dt:cursor_timevec(end))';
                cursor_pos_interp = zeros(length(tGrid),size(cursor_pos,2));
                for i=1:size(cursor_pos,2)
                    cursor_pos_interp(:,i) = interp1(cursor_timevec,cursor_pos(:,i),tGrid);
                end
                
                cursor_timevec = tGrid;
                cursor_pos = cursor_pos_interp;
            else
                dt = mean(diff(cursor_timevec));
            end

            % resample signals to match new time vector (assumes a uniform sampling rate)
            [temp_resampled,t_resamp] = resample_signals(cursor_pos,cursor_timevec,struct( ...
                'bin_size',bin_size, ...
                'samprate',1/dt,...
                'filter_n',filter_n,...
                'kaiser_beta',kaiser_beta...
                ));
            cursor_pos_resamp = interp1(t_resamp,temp_resampled,trial.timevec);

            % replace old cursor pos with this one where applicable
            error_cursor_idx = ~isnan(cursor_pos_resamp(:,1));
            trial.cursor_pos(error_cursor_idx,:) = cursor_pos_resamp(error_cursor_idx,:);
            trial.rel_cursor_pos(error_cursor_idx,:) = trial.cursor_pos(error_cursor_idx,:)-trial.ct_location(1:2);
            
            trial.cst_cursor_command(error_cursor_idx,1) = trial.lambda*(trial.rel_hand_pos(error_cursor_idx,1)+trial.rel_cursor_pos(error_cursor_idx,1));
            trial.cst_cursor_command(error_cursor_idx,2) = zeros(sum(error_cursor_idx),1);
            
            % get start and end indices (not the same as go and reward times)
            trial.idx_cstStartTime = find(error_cursor_idx,1,'first');
            trial.idx_cstEndTime = find(error_cursor_idx,1,'last');
        end
    else
        trial.hand_pos = [];
        trial.cursor_pos = [];
        trial.rel_hand_pos = [];
        trial.rel_cursor_pos = [];
        trial.cst_cursor_command = [];
    end
end

function trial = parse_smile_eye_data(in_trial,smile_trial,params)
    % parse out eye tracking (eye position and pupil diameter) from smile data struct
    % assumes in_trial has timevec based on the endTime found in parse_smile_events

    assert(isfield(in_trial,'timevec'),'Missing temporary timevector in trial data for some reason...')

    % assert length of analog channels divided by 1000 matches the end time of timevec
    assert(...
        size(smile_trial.TrialData.analogData,1)/double(1e3)==in_trial.timevec(end),...
        'Length of analog channel data does not match length of trial time...'...
        )

    % assume we're going for millisecond binning
    bin_size = 0.001;
    fill_err = false;
    assignParams(who,params)
    
    % copy over input trial
    trial = in_trial;

    % check to make sure eye tracking data is there
    eye_channel_names = {...
        'Left Eye X',...
        'Left Eye Y',...
        'Right Eye X',...
        'Right Eye Y',...
        'Left Pupil',...
        'Right Pupil'...
    };
    all_channel_names = smile_trial.Definitions.analogChannelNames;
    if ~all(ismember(eye_channel_names,all_channel_names))
        warning('Missing some or all of the eye data channels')
    end

    for eye_channel_num = 1:length(eye_channel_names)
        chan_idx = strcmpi(all_channel_names,eye_channel_names{eye_channel_num});
    end
end

function trial = parse_smile_spikes(in_trial,smile_trial,unit_guide,params)
    % parse out neural activity from smile data struct
    
    % import parameters
    array_alias = {};
    assignParams(who,params);

    % copy over input trial
    trial = in_trial;

    % Extract info from smile_trial
    array_name = smile_trial.TrialData.TDT.brainArea;

    % alias array names if we want
    if ~isempty(array_alias)
        for array_idx = 1:size(array_alias,1)
            if strcmp(array_name,array_alias{1,1})
                array_name = array_alias{1,2};
            end
        end
    end

    array_name = strrep(array_name,' ','');

    % construct a matrix of binned spikes and fill with spikes
    binned_spikes = zeros(length(trial.timevec),size(unit_guide,1));

    % fill matrix with histcounts
    spike_times = double(smile_trial.TrialData.TDT.snippetInfo');
    if ~isempty(spike_times)
        t_bin_edges = [trial.timevec trial.timevec(end)+mode(diff(trial.timevec))];
        for unit_idx = 1:size(unit_guide,1)
            unit_to_bin = unit_guide(unit_idx,:);
            filter_idx = spike_times(:,1)==unit_to_bin(1) & spike_times(:,2)==unit_to_bin(2);
            binned_spikes(:,unit_idx) = histcounts(spike_times(filter_idx,3)/1000,t_bin_edges)';
        end
    end

    % add to trial
    trial.(sprintf('%s_unit_guide',array_name)) = unit_guide;
    trial.(sprintf('%s_spikes',array_name)) = binned_spikes;
end

function unit_guide = get_smile_unit_guide(smile_data,params)
    % function to extract unit guide from all trials of smile_data
    
    % set up params
    valid_sort_idx = 1:10;
    assignParams(who,params);

    % concatenate snippet info from all trials of smile_data
    TrialData = horzcat(smile_data.TrialData);
    for trialnum = 1:length(TrialData)
        TrialData(trialnum).TDT.timestampUnits = 'ms';
    end
    TDT = horzcat(TrialData.TDT);
    snippetInfo = horzcat(TDT.snippetInfo);

    % get unit guide
    spike_times = double(snippetInfo');
    if ~isempty(spike_times)
        unit_guide = unique(spike_times(:,1:2),'rows');
        keep_units = ismember(unit_guide(:,2),valid_sort_idx);
        unit_guide = unit_guide(keep_units,:);
    else
        unit_guide = [-1 -1];
    end
end
