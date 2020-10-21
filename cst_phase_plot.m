%% Set up
    dataroot = '/data/raeed/project-data/smile/cst-gainlag';
    if ispc
        dataroot = 'C:\Users\Raeed\data\project-data\smile\cst-gainlag';
    end

    file_info = dir(fullfile(dataroot,'library','*COCST*'));
    filenames = horzcat({file_info.name})';
    
%% set up initial variables
    max_trial_length = 6;
    cursor_max = 50;
    
%% Loop through files to check out temporal smoothness at different lambdas
filetic = tic;
control_params_file_cell = cell(length(filenames),1);
for filenum = 31%1:length(filenames)
    td = load(fullfile(dataroot,'library',filenames{filenum}));
    td = td.trial_data;
    
%     % get rid of unsorted neurons
%     bad_units = td(1).M1_unit_guide(:,2)<=1;
%     % for file 31 only! remove neurons 8,2 and 64,2
%     if contains(td(1).date_time,'2018/06/26')
%         corr_units = [8 2;64 2];
%         bad_units = bad_units | ismember(td(1).M1_unit_guide,corr_units,'rows');
%     end
%     for trialnum = 1:length(td)
%         td(trialnum).M1_spikes = td(trialnum).M1_spikes(:,~bad_units);
%         td(trialnum).M1_unit_guide = td(trialnum).M1_unit_guide(~bad_units,:);
%     end
%     
%     if isempty(td(1).M1_unit_guide)
%         fprintf('Skipping file %d because no spike data...\n',filenum)
%         continue
%     end
%     
%     [~,td_cst] = getTDidx(td,'task','CST','result','R');
    [~,td_cst] = getTDidx(td,'task','CST');
    if ~isempty(td_cst)
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
        cutoff_freq = 70; %Hz
        samp_rate = 1/td_cst(1).bin_size;
        [filt_b,filt_a] = butter(16,cutoff_freq/(samp_rate/2));
        td_cst = filterSignals(td_cst,struct('signals','hand_pos','filt_a',filt_a,'filt_b',filt_b));
        td_cst = getDifferential(td_cst,struct('signals','hand_pos','alias','hand_vel'));
        td_cst = getDifferential(td_cst,struct('signals','hand_vel','alias','hand_acc'));
        td_cst = getDifferential(td_cst,struct('signals','cursor_pos','alias','cursor_vel'));
        
        % reorder for niceness
        td_cst = reorderTDfields(td_cst);
        
        % trim to 0.5s-5.5s 
%         td_cst = trimTD(td_cst,{'idx_goCueTime',floor(0.5/td_cst(1).bin_size)},{'idx_goCueTime',floor(5.5/td_cst(1).bin_size)});
        % trim to only CST window
        td_cst = trimTD(td_cst,'idx_cstStartTime','idx_cstEndTime');
        
        % get tolerable instability over time for each trial
        for trialnum = 1:length(td_cst)
            trial_length = size(td_cst(trialnum).cursor_pos,1);
            timevec = (0:(trial_length-1))'*td_cst(trialnum).bin_size;
            time_left = max_trial_length-timevec;
            td_cst(trialnum).tol_instab = 1./time_left .* log(cursor_max./abs(td_cst(trialnum).cursor_pos(:,1)+td_cst(trialnum).hand_pos(:,1)));
        end

        % plot out hand pos vs cursor pos for individual trials
        figure('defaultaxesfontsize',18)
        for trialnum = 1:length(td_cst)
            clf
            subplot(1,2,1)
            plot([-60 60],[0 0],'-k','linewidth',1)
            hold on
            plot([0 0],[-60 60],'-k','linewidth',1)
            plot([-60 60],[60 -60],'-k','linewidth',1)
            plot([-50 -50],[-60 60],'-r','linewidth',1)
            if strcmpi(td_cst(trialnum).result,'R')
                plot([-50 -50],[-60 60],'-g','linewidth',1)
                plot([50 50],[-60 60],'-g','linewidth',1)
            else
                plot([-50 -50],[-60 60],'-r','linewidth',1)
                plot([50 50],[-60 60],'-r','linewidth',1)
            end
            
%             scatter(td_cst(trialnum).cursor_pos(:,1),td_cst(trialnum).hand_pos(:,1),[],1:length(td_cst(trialnum).cursor_pos),'filled')
            scatter(td_cst(trialnum).cursor_pos(:,1),td_cst(trialnum).hand_pos(:,1),[],td_cst(trialnum).tol_instab,'filled')
            colorbar

%             plot3(td_cst(trialnum).cursor_pos(:,1),td_cst(trialnum).hand_pos(:,1),td_cst(trialnum).hand_vel(:,1),'-k')

            colormap(viridis)
            
            axis equal
            set(gca,'box','off','tickdir','out','xlim',[-60 60],'ylim',[-60 60],'clim',[0 5])
            title(strcat(sprintf('%s %s \\lambda = %f, Trial ID: %d',...
                td_cst(trialnum).monkey,...
                td_cst(trialnum).date_time,...
                td_cst(trialnum).lambda,...
                td_cst(trialnum).trial_id)))
            xlabel('Cursor position (mm)')
            ylabel('Hand position (mm)')
            
            subplot(1,2,2)
            plot([0 6],td_cst(trialnum).lambda*[1 1],'k','linewidth',1)
            hold on
            scatter(td_cst(trialnum).bin_size*(1:length(td_cst(trialnum).cursor_pos)),td_cst(trialnum).tol_instab,[],td_cst(trialnum).tol_instab,'filled')
            set(gca,'box','off','tickdir','out','ylim',[0 5],'clim',[0 5])
            waitforbuttonpress;
        end
        
        fprintf('Finished file %d of %d at time %f\n',filenum,length(filenames),toc(filetic))
    else
        fprintf('Incomplete dataset for file %d\n',filenum)
    end
    control_params = vertcat(control_params_file_cell{:});
end