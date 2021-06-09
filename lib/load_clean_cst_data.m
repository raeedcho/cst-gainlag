function [td_cst] = load_clean_cst_data(filename,params)
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
%       td_cst - TrialData struture with clean CST trials in it

    % set up params
    cutoff_freq = 70; %Hz <- this is the cutoff to get rid of the weird harmonics in velocity and acceleration
    if nargin>1
        assignParams(who,params);
    end
    
    % load file
    td = load(filename);
    td = td.trial_data;
    [~,td_cst] = getTDidx(td,'task','CST');
    if isempty(td_cst)
        return
    end

    % remove aborts
    abort_idx = isnan(cat(1,td_cst.idx_goCueTime));
    td_cst(abort_idx) = [];
    fprintf('Removed %d trials that monkey aborted\n',sum(abort_idx))

    % trim nans off the end of trials for cursor pos and hand pos
    for trialnum = 1:length(td_cst)
        nan_times = any(isnan(getSig(td_cst(trialnum),getTDfields(td_cst,'time'))),2);
        first_viable_time = find(~nan_times,1,'first');
        last_viable_time = find(~nan_times,1,'last');
        td_cst(trialnum) = trimTD(td_cst(trialnum),{'start',first_viable_time-1},{'start',last_viable_time-1});
    end

    % fill in CST windows (go cue and reward time seem to be off by some random amount)
    for trialnum = 1:length(td_cst)
        cst_window = td_cst(trialnum).cursor_pos(:,1)~=td_cst(trialnum).hand_pos(:,1);
        td_cst(trialnum).idx_cstStartTime = find(cst_window,1,'first');
        td_cst(trialnum).idx_cstEndTime = find(cst_window,1,'last');
    end

    % fill kinematic signals (filter ahead of differentials)
    samp_rate = 1/td_cst(1).bin_size;
    [filt_b,filt_a] = butter(16,cutoff_freq/(samp_rate/2));
    td_cst = filterSignals(td_cst,struct('signals','hand_pos','filt_a',filt_a,'filt_b',filt_b));
    td_cst = getDifferential(td_cst,struct('signals','hand_pos','alias','hand_vel'));
    td_cst = getDifferential(td_cst,struct('signals','hand_vel','alias','hand_acc'));
    td_cst = getDifferential(td_cst,struct('signals','cursor_pos','alias','cursor_vel'));

    % add shifted signals
    shift_bins = floor(-0.1/td_cst(1).bin_size);
    td_cst = dupeAndShift(td_cst,'cursor_pos',shift_bins,'cursor_vel',shift_bins);

    % reorder for niceness
    td_cst = reorderTDfields(td_cst);
