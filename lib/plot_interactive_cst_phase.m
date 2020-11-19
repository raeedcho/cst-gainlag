function [handle] = plot_interactive_cst_phase(td_cst)
% plot out hand pos vs cursor pos for individual trials
% use 'h' and 'l' keys to go back and forward, respectively
% 'q' to quit

    handle = figure('defaultaxesfontsize',18,'position',[1544 397 813 839]);
    colormap(viridis)
    trialnum = 1;
    while trialnum<=length(td_cst)
        clf
        phase_ax = subplot(4,4,[9 10 13 14]);
        plot_cst_phase(td_cst(trialnum))
        xlabel('Cursor position (mm)')
        ylabel('Hand position (mm)')

        hand_ax = subplot(4,4,[11 12 15 16]);
        plot_cst_signal(td_cst(trialnum),struct('sig_to_plot','hand_pos'));
        set(gca,'ylim',60*[-1 1])
        xlabel('Time after CST start (s)')

        cursor_ax = subplot(4,4,[1 2 5 6]);
        plot_cst_signal(td_cst(trialnum),struct('sig_to_plot','cursor_pos','flipxy',true))
        if strcmpi(td_cst(trialnum).result,'R')
            plot([-50 -50],[0 6],'-g','linewidth',1)
            plot([50 50],[0 6],'-g','linewidth',1)
        else
            plot([-50 -50],[0 6],'-r','linewidth',1)
            plot([50 50],[0 6],'-r','linewidth',1)
        end
        set(gca,'xlim',60*[-1 1])
        ylabel('Time after CST start (s)')

        hand_vel_ax = subplot(4,4,[7 8]);
%         plot_cst_signal(td_cst(trialnum),struct('sig_to_plot',{{'hand_pos',2}}))
        scatter(...
            td_cst(trialnum).bin_size*(1:length(td_cst(trialnum).cursor_pos)),...
            -td_cst(trialnum).hand_pos(:,2),...
            [],1:length(td_cst(trialnum).cursor_pos),'filled')
        set(gca,'box','off','tickdir','out')

        hand_acc_ax = subplot(4,4,[3 4]);
%         plot_cst_signal(td_cst(trialnum),struct('sig_to_plot','hand_acc'))
        for restorenum = 1:length(td_cst(trialnum).idx_restoreStart)
            blockstart = td_cst(trialnum).idx_restoreStart(restorenum);
            blockend = td_cst(trialnum).idx_restoreEnd(restorenum);
            patch(...
                td_cst(trialnum).bin_size*[blockstart blockend blockend blockstart],...
                60*[-1 -1 1 1],...
                [0.8 0.8 0.8],...
                'edgecolor','none')
            hold on
        end
        plot([0 6],[0 0],'k','linewidth',1)
        scatter(...
            td_cst(trialnum).bin_size*(1:length(td_cst(trialnum).cursor_pos)),...
            (td_cst(trialnum).hand_pos(:,1)+td_cst(trialnum).cursor_pos(:,1)),...
            [],1:length(td_cst(trialnum).cursor_pos),'filled')
        set(gca,'box','off','tickdir','out')

        % link axes for interactivity
        linkaxes([phase_ax cursor_ax],'x')
        linkaxes([phase_ax hand_ax],'y')
        linkaxes([hand_ax hand_vel_ax hand_acc_ax],'x')

        suptitle(strcat(sprintf('%s %s \\lambda = %f, Trial ID: %d',...
            td_cst(trialnum).monkey,...
            td_cst(trialnum).date_time,...
            td_cst(trialnum).lambda,...
            td_cst(trialnum).trial_id)))

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
function plot_cst_phase(trial,params)
% plot hand position against cursor position, plus all the extras

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
        trial.cursor_pos(:,1),...
        trial.hand_pos(:,1),...
        [],1:length(trial.cursor_pos),'filled')

    axis equal
    set(gca,'box','off','tickdir','out','xlim',[-60 60],'ylim',[-60 60])
end

function plot_cst_signal(trial,params)
% plot a signal of the CST trial

    flipxy = false;
    sig_to_plot = '';
    assignParams(who,params)
    
    sig_data = getSig(trial,sig_to_plot);
    time_data = trial.bin_size*(1:length(trial.cursor_pos));
    
    if flipxy
        % plot restoration blocks
        for restorenum = 1:length(trial.idx_restoreStart)
            blockstart = trial.idx_restoreStart(restorenum);
            blockend = trial.idx_restoreEnd(restorenum);
            patch(...
                60*[-1 -1 1 1],...
                trial.bin_size*[blockstart blockend blockend blockstart],...
                [0.8 0.8 0.8],...
                'edgecolor','none')
            hold on
        end
        plot([0 0],[0 6],'k','linewidth',1)
        scatter(...
            sig_data(:,1),...
            time_data,...
            [],1:length(trial.cursor_pos),'filled')
    else
        % plot restoration blocks
        for restorenum = 1:length(trial.idx_restoreStart)
            blockstart = trial.idx_restoreStart(restorenum);
            blockend = trial.idx_restoreEnd(restorenum);
            patch(...
                trial.bin_size*[blockstart blockend blockend blockstart],...
                60*[-1 -1 1 1],...
                [0.8 0.8 0.8],...
                'edgecolor','none')
            hold on
        end
        plot([0 6],[0 0],'k','linewidth',1)
        scatter(...
            time_data,...
            sig_data(:,1),...
            [],1:length(trial.cursor_pos),'filled')
    end
    set(gca,'box','off','tickdir','out')
end
    
    