function td_cst = calcTolInstab(td_cst,params)
% calculate time-varying tolerable instability in trial for CST task
% Adds 'tol_instab' signal to trial data structure in output
%
% Inputs:
%   td_cst - TrialData structure with only CST trials in it
%   params - list of parameters for running
%       .max_trial_length - maximum length of trial in seconds (default: 6)
%       .cursor_max - maximum bounds of cursor position in mm (default: 50)
%       .assume_zero_hand_pos - set true if we want to assume hand pos is zero;
%           set false if we want to assume hand pos stays still at current position
%           (default: false)
%
% Outputs:
%   td_cst - TrialData structure with new 'tol_instab' signal added to it

    % params
    max_trial_length = 6;
    cursor_max = 50;
    assume_zero_hand_pos = false;
    if nargin>1
        assignParams(who,params);
    end

    % get tolerable instability over time for each trial
    for trialnum = 1:length(td_cst)
        trial_length = size(td_cst(trialnum).cursor_pos,1);
        timevec = (0:(trial_length-1))'*td_cst(trialnum).bin_size;
        time_left = max_trial_length-timevec;
        if assume_zero_hand_pos
            td_cst(trialnum).tol_instab = 1./time_left .* log(cursor_max./abs(td_cst(trialnum).cursor_pos(:,1)));
        else
            td_cst(trialnum).tol_instab = 1./time_left .* log(cursor_max./abs(td_cst(trialnum).cursor_pos(:,1)+td_cst(trialnum).hand_pos(:,1)));
        end
    end