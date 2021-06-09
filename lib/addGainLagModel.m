function td = addGainLagModel(td,params)
% adds a control signal proportional to a lagged cursor position

    cursor_sig = {'cursor_sig',1};
    out_sig_name = 'gainlag_model';
    gain = 1;
    lag = 0;
    if nargin>1
        assignParams(who,params)
    end
    
    % check signals
    cursor_sig = check_signals(td,cursor_sig);
    
    % check lag
    assert(lag>=0,'Currently, lag must be >=0')
    lag_bins = floor(lag/td(1).bin_size);
    if lag_bins*td(1).bin_size~=lag
        warning('Lag doesn''t quite match bin size. Actual lag will be %f seconds\n',lag_bins*td(1).bin_size)
    end
    
    % go through trials
    for trialnum = 1:length(td)
        cursor_data = get_vars(td(trialnum),cursor_sig);
        out_sig = nan(size(cursor_data));
        
        out_sig(lag_bins+1:end,:) = gain*cursor_data(1:end-lag_bins,:);
        
        td(trialnum).(out_sig_name) = out_sig;
    end