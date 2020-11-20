function [handle] = plot_interactive_cst_phase(td_cst,params)
% plot out hand pos vs cursor pos for individual trials
% use 'h' and 'l' keys to go back and forward, respectively
% 'q' to quit

    % set up figure to save time later
    handle = figure(...
        'defaultaxesfontsize',18,...
        'position',[1544 397 813 839],...
        'NextPlot','add',...
        'Name',sprintf('%s %s',td_cst(1).monkey,td_cst(1).date_time));
%     colormap(viridis)
    colormap([1 1 1;0.8 0.8 0.8])
%     colormap(flipud(gray))
    phase_ax = subplot(4,4,[9 10 13 14]);
    hand_ax = subplot(4,4,[11 12 15 16]);
    cursor_ax = subplot(4,4,[1 2 5 6]);
    extra1_ax = subplot(4,4,[3 4]);
    extra2_ax = subplot(4,4,[7 8]);
    
    set(phase_ax,'box','off','tickdir','out','nextplot','add','xlim',[-60 60],'ylim',[-60 60],'dataaspectratio',[1 1 1])
    set(hand_ax,'box','off','tickdir','out','nextplot','add','xlim',[0 6],'ylim',[-60 60])
    set(cursor_ax,'box','off','tickdir','out','nextplot','add','xlim',[-60 60],'ylim',[0 6])
    set(extra1_ax,'box','off','tickdir','out','nextplot','add','xlim',[0 6])
    set(extra2_ax,'box','off','tickdir','out','nextplot','add','xlim',[0 6])
    
    linkaxes([phase_ax cursor_ax],'x')
    linkaxes([phase_ax hand_ax],'y')
    linkaxes([hand_ax extra1_ax extra2_ax],'x')
    
    trialnum = 1;
    while trialnum<=length(td_cst)
        % clear axes
        cla(phase_ax)
        cla(hand_ax)
        cla(cursor_ax)
        cla(extra1_ax)
        cla(extra2_ax)
        
        plot_cst_phase(phase_ax,td_cst(trialnum))
        xlabel(phase_ax,'Cursor position (mm)')
        ylabel(phase_ax,'Hand position (mm)')

        plot_cst_signal(hand_ax,td_cst(trialnum),struct('sig_to_plot','hand_pos'));
        xlabel(hand_ax,'Time after CST start (s)')

        plot_cst_signal(cursor_ax,td_cst(trialnum),struct('sig_to_plot','cursor_pos','flipxy',true))
        if strcmpi(td_cst(trialnum).result,'R')
            plot(cursor_ax,[-50 -50],[0 6],'-g','linewidth',1)
            plot(cursor_ax,[50 50],[0 6],'-g','linewidth',1)
        else
            plot(cursor_ax,[-50 -50],[0 6],'-r','linewidth',1)
            plot(cursor_ax,[50 50],[0 6],'-r','linewidth',1)
        end
        ylabel(cursor_ax,'Time after CST start (s)')

%         plot_cst_signal(extra1_ax,td_cst(trialnum),struct('sig_to_plot',{{'hand_pos',2}}))
        scatter(extra1_ax,...
            td_cst(trialnum).bin_size*(1:length(td_cst(trialnum).cursor_pos)),...
            -td_cst(trialnum).hand_pos(:,2),...
            [],viridis(length(td_cst(trialnum).cursor_pos)),'filled')

%         plot_cst_signal(td_cst(trialnum),struct('sig_to_plot','hand_acc'))
        for restorenum = 1:length(td_cst(trialnum).idx_restoreStart)
            blockstart = td_cst(trialnum).idx_restoreStart(restorenum);
            blockend = td_cst(trialnum).idx_restoreEnd(restorenum);
            patch(extra2_ax,...
                td_cst(trialnum).bin_size*[blockstart blockend blockend blockstart],...
                60*[-1 -1 1 1],...
                [0.8 0.8 0.8],...
                'edgecolor','none')
            hold on
        end
        plot(extra2_ax,[0 6],[0 0],'k','linewidth',1)
        scatter(extra2_ax,...
            td_cst(trialnum).bin_size*(1:length(td_cst(trialnum).cursor_pos)),...
            (td_cst(trialnum).hand_pos(:,1)+td_cst(trialnum).cursor_pos(:,1)),...
            [],viridis(length(td_cst(trialnum).cursor_pos)),'filled')

        % link axes for interactivity (do it manually to save time)
%         drawnow % flush buffer to make sure we have the right lims
%         linkaxes([phase_ax cursor_ax],'x')
%         linkaxes([phase_ax hand_ax],'y')
%         linkaxes([hand_ax hand_vel_ax hand_acc_ax],'x')
% 
%         title(cursor_ax,strcat(sprintf('\\lambda = %f, Trial ID: %d',...
%             td_cst(trialnum).lambda,...
%             td_cst(trialnum).trial_id)))

        % set up navigation keys
        while true
            if waitforbuttonpress==1
                charpressed = get(gcf,'CurrentCharacter');
                if charpressed == 'h'
                    trialnum = max(trialnum-1,1);
                    break;
                elseif charpressed == 'l'
                    trialnum = min(trialnum+1,length(td_cst));
                    break;
                elseif charpressed == 'q'
                    return
                end
            end
        end
    end
end

%%%%% Subfunctions
function plot_cst_phase(ax,trial)
% plot hand position against cursor position, plus all the extras

    plot(ax,[-60 60],[0 0],'-k','linewidth',1)
    hold on
    plot(ax,[0 0],[-60 60],'-k','linewidth',1)
    plot(ax,[-60 60],[60 -60],'-k','linewidth',1)
    plot(ax,[-50 -50],[-60 60],'-r','linewidth',1)
    patch(ax,...
        [0 0 -60 60 0],...
        [-60 60 60 -60 -60],...
        [0.8 0.8 0.8],'edgecolor','none')
    if strcmpi(trial.result,'R')
        plot(ax,[-50 -50],[-60 60],'-g','linewidth',1)
        plot(ax,[50 50],[-60 60],'-g','linewidth',1)
    else
        plot(ax,[-50 -50],[-60 60],'-r','linewidth',1)
        plot(ax,[50 50],[-60 60],'-r','linewidth',1)
    end

%     scatter(ax,...
%         trial.cursor_pos(:,1),...
%         trial.hand_pos(:,1),...
%         [],viridis(length(trial.cursor_pos)),'filled')
    
    patch(ax,...
        [trial.cursor_pos(:,1);NaN],...
        [trial.hand_pos(:,1);NaN],...
        permute([viridis(length(trial.cursor_pos));nan(1,3)],[1 3 2]),...
        'facecolor','none',...
        'edgecolor','interp',...
        'marker','.',...
        'linewidth',5)
    
end

function plot_cst_signal(ax,trial,params)
% plot a signal of the CST trial

    flipxy = false;
    sig_to_plot = '';
    patch_sig = 'is_restoring';
    assignParams(who,params)
    
    sig_data = getSig(trial,sig_to_plot);
    time_data = trial.bin_size*(1:length(trial.cursor_pos))';
    patch_data = getSig(trial,patch_sig);
    
    %set up patch
    vertices = horzcat(repmat(time_data,2,1),60*vertcat(-ones(length(time_data),1),ones(length(time_data),1)));
    faces = [1:length(time_data) (length(time_data)*2):-1:(length(time_data)+1)];
    vert_colors = repmat(patch_data,2,1);
    
    if flipxy
        patch(ax,...
            'faces',faces,...
            'vertices',vertices(:,[2 1]),...
            'facevertexcdata',vert_colors,...
            'facecolor','interp',...
            'edgecolor','none')
        hold on
        plot(ax,[0 0],[0 6],'k','linewidth',1)
%         scatter(ax,...
%             sig_data(:,1),...
%             time_data,...
%             [],viridis(length(trial.cursor_pos)),'filled')
        patch(ax,...
            [sig_data(:,1);NaN],...
            [time_data;NaN],...
            permute([viridis(length(sig_data));nan(1,3)],[1 3 2]),...
            'facecolor','none',...
            'edgecolor','interp',...
            'linewidth',5)
    else
        patch(ax,...
            'faces',faces,...
            'vertices',vertices,...
            'facevertexcdata',vert_colors,...
            'facecolor','interp',...
            'edgecolor','none')
        hold on
        
        plot(ax,[0 6],[0 0],'k','linewidth',1)
%         scatter(ax,...
%             time_data,...
%             sig_data(:,1),...
%             [],viridis(length(trial.cursor_pos)),'filled')
        patch(ax,...
            [time_data;NaN],...
            [sig_data(:,1);NaN],...
            permute([viridis(length(sig_data));nan(1,3)],[1 3 2]),...
            'facecolor','none',...
            'edgecolor','interp',...
            'linewidth',5)
    end
end
    
    