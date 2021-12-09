function plot_cst_traces(trial,params)
% plot hand and cursor traces of the CST trial

    cursor_sig = {'rel_cursor_pos',1};
    hand_sig = {'rel_hand_pos',1};
    assignParams(who,params)
    
    time_data = trial.bin_size*(1:length(trial.cursor_pos));

    plot([0 6],[0 0],'k','linewidth',1)
    hold on
    plot(...
        time_data,...
        getSig(trial,cursor_sig),...
        '-b','linewidth',2)
    plot(...
        time_data,...
        getSig(trial,hand_sig),...
        '-r','linewidth',2)
    set(gca,'box','off','tickdir','out')
end
    