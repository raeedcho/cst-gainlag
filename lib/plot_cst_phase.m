function plot_cst_phase(trial,params)
% plot hand position against cursor position, plus all the extras

    cursor_sig = {'cursor_pos',1};
    hand_sig = {'hand_pos',1};
    color_sig = {};
    assignParams(who,params)
    
    if isempty(color_sig)
        color_vals = 1:length(get_vars(trial,cursor_sig));
    else
        color_vals = get_vars(trial,color_sig);
    end

    plot([-60 60],[0 0],'-k','linewidth',1)
    hold on
    plot([0 0],[-60 60],'-k','linewidth',1)
    plot([-60 60],[60 -60],'-k','linewidth',1)
    plot([-50 -50],[-60 60],'-r','linewidth',1)
    patch(...
        [0 0 -60 60 0],...
        [-60 60 60 -60 -60],...
        [0.8 0.8 0.8],'edgecolor','none')
    if strcmpi(trial.result,'R')
        plot([-50 -50],[-60 60],'-g','linewidth',1)
        plot([50 50],[-60 60],'-g','linewidth',1)
    else
        plot([-50 -50],[-60 60],'-r','linewidth',1)
        plot([50 50],[-60 60],'-r','linewidth',1)
    end

    scatter(...
        get_vars(trial,cursor_sig),...
        get_vars(trial,hand_sig),...
        [],color_vals,'filled')

    axis equal
    set(gca,'box','off','tickdir','out','xlim',[-60 60],'ylim',[-60 60])
end