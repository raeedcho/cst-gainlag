%% Set up
    dataroot = '/data/raeed/project-data/smile/cst-gainlag';
    if ispc
        dataroot = 'C:\Users\Raeed\data\project-data\smile\cst-gainlag';
    end

    file_info = dir(fullfile(dataroot,'library','*COCST*'));
    filenames = horzcat({file_info.name})';

%% Loop through files to check out temporal smoothness at different lambdas
filetic = tic;
control_params_file_cell = cell(length(filenames),1);
for filenum = 1:length(filenames)
    td = load(fullfile(dataroot,'library',filenames{filenum}));
    td = td.trial_data;
    
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
    
    if isempty(td(1).M1_unit_guide)
        fprintf('Skipping file %d because no spike data...\n',filenum)
        continue
    end
    
    [~,td_cst] = getTDidx(td,'task','CST','result','R');
    if ~isempty(td_cst)
        % remove aborts
        abort_idx = isnan(cat(1,td_cst.idx_goCueTime));
        td_cst(abort_idx) = [];
        fprintf('Removed %d trials that monkey aborted\n',sum(abort_idx))
        
        % early trim to remove nans on end for filtering
        td_cst = trimTD(td_cst,{'idx_startTime',0},{'idx_rewardTime',floor(0.13/td_cst(1).bin_size)});
        
        % fill kinematic signals (filter ahead of differentials)
        cutoff_freq = 70; %Hz
        samp_rate = 1/td_cst(1).bin_size;
        [filt_b,filt_a] = butter(16,cutoff_freq/(samp_rate/2));
        td_cst = filterSignals(td_cst,struct('signals','hand_pos','filt_a',filt_a,'filt_b',filt_b));
        td_cst = getDifferential(td_cst,struct('signals','hand_pos','alias','hand_vel'));
        td_cst = getDifferential(td_cst,struct('signals','hand_vel','alias','hand_acc'));
        
        % trim to 0.5s-5.5s 
        td_cst = trimTD(td_cst,{'idx_goCueTime',floor(0.5/td_cst(1).bin_size)},{'idx_goCueTime',floor(5.5/td_cst(1).bin_size)});

        td_cst = removeBadNeurons(td_cst,struct(...
            'arrays',{{'M1'}},...
            'min_fr',0.1,...
            'calc_fr',true,...
            'do_shunt_check',true));
        
        % collect trial control information
        control_params_cell = cell(length(td_cst),1);
        for trialnum = 1:length(td_cst)
            % Things to collect:
            % - monkey
            % - date
            % - trial ID
            % - trial result
            % - lambda
            % - cursor start
            % - trial length (?)
            % - fitted control parameters
            % - dynamics eigenvalues
            % - max positive real eigenvalue
            % - dynamics eigenvectors
            % - eigenvalue bound
            
            % calculate control parameters
            hand_cursor_state = [td_cst(trialnum).cursor_pos(:,1) td_cst(trialnum).hand_pos(:,1) td_cst(trialnum).hand_vel(:,1)];
            control_k = hand_cursor_state\td_cst(trialnum).hand_acc(:,1);
            control_k = control_k';
            
            % Calculate controlled eigenvalues
            A = [...
                td_cst(trialnum).lambda td_cst(trialnum).lambda 0;...
                0 0 1;...
                0 0 0];
            B = [0;0;1];
            [eigvec,eigval] = eig(A+B*control_k);
            
            % calculate maximum possible positive eigenvalue given start position
            max_trial_length = 5.5;
            cursor_max = 50;
            cursor_start = td_cst(trialnum).cursor_pos(1,1);
            eigval_bound = 1/max_trial_length * log(cursor_max/abs(cursor_start));
            
            control_params_cell{trialnum} = table(...
                {td_cst(trialnum).monkey},...
                {td_cst(trialnum).date_time},...
                td_cst(trialnum).trial_id,...
                {td_cst(trialnum).result},...
                td_cst(trialnum).lambda,...
                cursor_start,...
                size(hand_cursor_state,1)*td_cst(trialnum).bin_size,...
                control_k,...
                diag(eigval)',...
                max(real(diag(eigval))),...
                {eigvec},...
                eigval_bound,...
                'VariableNames',{...
                    'monkey',...
                    'date_time',...
                    'trial_id',...
                    'result',...
                    'lambda',...
                    'cursor_start',...
                    'trial_length',...
                    'control_k',...
                    'eigval',...
                    'max_real_eigval',...
                    'eigvec',...
                    'eigval_bound'});
            control_params_cell{trialnum}.Properties.VariableDescriptions = {'meta','meta','meta','meta','meta','linear','linear','linear','linear','linear','linear','linear'};
        end
        control_params_file_cell{filenum} = vertcat(control_params_cell{:});
        
        weird_trials = find(control_params_file_cell{filenum}.eigval_bound < control_params_file_cell{filenum}.max_real_eigval);
        fprintf('There are %d weird trials out of %d in file %d\n',length(weird_trials),length(td_cst),filenum)
        
        % plot out hand pos vs cursor pos for individual trials
%         figure
%         for trialnum_iter = 1:length(weird_trials)
%             trialnum = weird_trials(trialnum_iter);
% %         for trialnum = 1:length(td_cst)
%             clf
% %             plot([-60 60],[0 0],'-k','linewidth',1)
% %             hold on
% %             plot([0 0],[-60 60],'-k','linewidth',1)
% %             plot([-60 60],[60 -60],'-k','linewidth',1)
% %             plot(td_cst(trialnum).cursor_pos(:,1),td_cst(trialnum).hand_pos(:,1),'-ok','linewidth',2)
%             plot3(td_cst(trialnum).cursor_pos(:,1),td_cst(trialnum).hand_pos(:,1),td_cst(trialnum).hand_vel(:,1),'-k')
% %             axis equal
% %             set(gca,'box','off','tickdir','out','xlim',[-60 60],'ylim',[-60 60])
%             title(strcat('\lambda = ',num2str(td_cst(trialnum).lambda),sprintf(', Trial id: %d', td_cst(trialnum).trial_id)))
%             xlabel('Cursor position (mm)')
%             ylabel('Hand position (mm)')
%             waitforbuttonpress;
%         end
        
        fprintf('Finished file %d of %d at time %f\n',filenum,length(filenames),toc(filetic))
    else
        fprintf('Incomplete dataset for file %d\n',filenum)
    end
    control_params = vertcat(control_params_file_cell{:});
end

%% Select an intermediate lambda and plot progress
lambda_specific = 2.4;
[~,control_specific] = getNTidx(control_params,'monkey','Ford','lambda',lambda_specific);

% gainlag_specific(gainlag_specific.gain<1,:) = [];
% they should already be sorted, so just plot out eigenvalues
figure('defaultaxesfontsize',18)
% subplot(2,1,1)
scatter(1:height(control_specific),control_specific.max_real_eigval,[],'k','filled')
xlabel('Trial')
ylabel('Max real eigenvalue')
set(gca,'box','off','tickdir','out')
% subplot(2,1,2)
% scatter(1:height(control_specific),control_specific.gain,[],'k','filled')
% xlabel('Trial')
% ylabel('Gain')
% set(gca,'box','off','tickdir','out')
title(sprintf('lambda = %0.1f',lambda_specific))

%% all trials
figure('defaultaxesfontsize',18)
subplot(2,1,1)
scatter(1:height(control_params),control_params.max_real_eigval,[],control_params.lambda,'filled')
colormap(viridis)
c = colorbar;
c.Label.String = '\lambda';
xlabel('Trial number')
ylabel('Max real eigenvalue')
set(gca,'box','off','tickdir','out')

% 3d?
figure('defaultaxesfontsize',18)
scatter3(datenum(control_params.date_time),control_params.lambda,control_params.max_real_eigval,[],control_params.trial_id,'filled')
colormap(viridis)
c = colorbar;
c.Label.String = 'Trial number';
xlabel('Date number')
ylabel('\lambda')
zlabel('Max real eigenvalue')
set(gca,'box','off','tickdir','out')