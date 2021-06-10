function [td] = load_clean_cst_data(filename,params)
% loads cst data from filename and cleans it up. Cleanup includes:
%   - removing aborted trials
%   - trimming trials to get rid of NaNs in cursor pos and hand pos signals
%   - determining the window of CST based on when cursor and hand are different
%   - filter and differentiate kinematic signals
%
%   Inputs:
%       filename - name of file to load
%       params
%           .cutoff_freq - cutoff frequency for hand_pos signal filtering
%
%   Outputs:
%       td - TrialData struture with clean CO and CST trials in it

    % set up params
    cutoff_freq = 70; %Hz <- this is the cutoff to get rid of the weird harmonics in velocity and acceleration
    if nargin>1
        assignParams(who,params);
    end
    
    % load file
    td = load(filename);
    td = td.trial_data;
    
    % remove aborts
    abort_idx = isnan(cat(1,td.idx_goCueTime));
    td(abort_idx) = [];
    fprintf('Removed %d trials that monkey aborted\n',sum(abort_idx))

    % trim nans off the end of trials for cursor pos and hand pos
    for trialnum = 1:length(td)
        nan_times = any(isnan(getSig(td(trialnum),getTDfields(td,'time'))),2);
        first_viable_time = find(~nan_times,1,'first');
        last_viable_time = find(~nan_times,1,'last');
        td(trialnum) = trimTD(td(trialnum),{'start',first_viable_time-1},{'start',last_viable_time-1});
    end

    % fill in CST windows (go cue and reward time seem to be off by some random amount)
    cst_idx = getTDidx(td,'task','CST');
    for trialnum = 1:length(cst_idx)
        trial_idx = cst_idx(trialnum);
        cst_window = td(trial_idx).cursor_pos(:,1)~=td(trial_idx).hand_pos(:,1);
        td(trial_idx).idx_cstStartTime = find(cst_window,1,'first');
        td(trial_idx).idx_cstEndTime = find(cst_window,1,'last');
    end
    
    % fill in CO target directions
    for trialnum = 1:length(td)
        if td(trialnum).task=='CO'
            td(trial_idx).tgtDir = atan2d(...
                td_co(trialnum).ot_location(2)-td_co(trialnum).ct_location(2),...
                td_co(trialnum).ot_location(1)-td_co(trialnum).ct_location(1));
        else
            td(trialnum).tgtDir = NaN;
        end
    end
    
    % fill kinematic signals (filter ahead of differentials)
    samp_rate = 1/td(1).bin_size;
    [filt_b,filt_a] = butter(16,cutoff_freq/(samp_rate/2));
    td = filterSignals(td,struct('signals','hand_pos','filt_a',filt_a,'filt_b',filt_b));
    td = getDifferential(td,struct('signals','hand_pos','alias','hand_vel'));
    td = getDifferential(td,struct('signals','hand_vel','alias','hand_acc'));
    td = getDifferential(td,struct('signals','cursor_pos','alias','cursor_vel'));

    % add shifted signals
    shift_bins = floor(-0.1/td(1).bin_size);
    td = dupeAndShift(td,'cursor_pos',shift_bins,'cursor_vel',shift_bins);
    
    % neural data stuff
    if ~isempty(td(1).M1_unit_guide)
        % get rid of unsorted neurons
        bad_units = td(1).M1_unit_guide(:,2)<=1;
        % for file 31 only! remove neurons 8,2 and 64,2
        if contains(td(1).date_time,'2018/06/26')
            corr_units = [8 2;64 2];
            bad_units = bad_units | ismember(td(1).M1_unit_guide,corr_units,'rows');
        end
        for trialnum = 1:length(td)
            td(trialnum).M1_spikes = td(trialnum).M1_spikes(:,~bad_units);
            td(trialnum).M1_unit_guide = td(trialnum).M1_unit_guide(~bad_units,:);
        end

        td = removeBadNeurons(td,struct(...
            'do_fr_check',true,...
            'min_fr',0.1,...
            'fr_window',{{'idx_cstStartTime',0;'idx_cstEndTime',0}},...
            'calc_fr',true,...
            'use_trials',getTDidx(td,'task','CST')));
    end

    % reorder for niceness
    td = reorderTDfields(td);
