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
for filenum = 1:length(filenames)
    td_cst = load_clean_cst_data(fullfile(dataroot,'library',filenames{filenum}));
    
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
    if isempty(td_cst)
        fprintf('Incomplete dataset for file %d\n',filenum)
        continue
    end
    td_cst = trimTD(td_cst,'idx_cstStartTime','idx_cstEndTime');
    td_cst = calcTolInstab(td_cst);
%     td_cst = findRestorationBlocks(td_cst);

    % trim all trials to the same length
    td_cst = trimTD(td_cst,struct(...
        'idx_start',{{'idx_cstStartTime',0}},...
        'idx_end',{{'idx_cstStartTime',floor(3/td_cst(1).bin_size)}},...
        'remove_short',true));

%         td_cst = removeBadNeurons(td_cst,struct(...
%             'arrays',{{'M1'}},...
%             'min_fr',0.1,...
%             'calc_fr',true,...
%             'do_shunt_check',true));

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
        % - gain/lag parameters
        % - fitted control parameters
        % - dynamics eigenvalues
        % - max positive real eigenvalue
        % - dynamics eigenvectors
        % - eigenvalue bound

        % gain/lag analysis
        [corr_trace,corr_lag] = xcorr(td_cst(trialnum).hand_vel(:,1),td_cst(trialnum).cursor_vel(:,1));
        [~,lag_idx] = min(corr_trace);
        y_motion = sqrt(mean(td_cst(trialnum).hand_pos(:,2).^2));
        x_motion = sqrt(mean(td_cst(trialnum).hand_pos(:,1).^2));
        gain = sqrt(mean(td_cst(trialnum).hand_vel(:,1).^2))/...
            sqrt(mean(td_cst(trialnum).cursor_vel(:,1).^2));

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
            gain,...
            corr_lag(lag_idx)*td_cst(trialnum).bin_size,...
            x_motion,...
            y_motion,...
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
                'gain',...
                'lag',...
                'x_motion',...
                'y_motion',...
                'control_k',...
                'eigval',...
                'max_real_eigval',...
                'eigvec',...
                'eigval_bound'});
        control_params_cell{trialnum}.Properties.VariableDescriptions = {...
            'meta','meta','meta','meta','meta','linear','linear','linear','linear','linear','linear','linear','linear','linear','linear','linear'};
    end
    control_params_file_cell{filenum} = vertcat(control_params_cell{:});

    weird_trials = find(control_params_file_cell{filenum}.eigval_bound < control_params_file_cell{filenum}.max_real_eigval);
    fprintf('There are %d weird trials out of %d in file %d\n',length(weird_trials),length(td_cst),filenum)

    fprintf('Finished file %d of %d at time %f\n',filenum,length(filenames),toc(filetic))
        
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