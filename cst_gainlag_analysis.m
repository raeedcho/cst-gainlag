%% Set up
    dataroot = '/data/raeed/project-data/smile/cst-gainlag';
    if ispc
        dataroot = 'C:\Users\Raeed\data\project-data\smile\cst-gainlag';
    end

    file_info = dir(fullfile(dataroot,'library','*COCST*'));
    filenames = horzcat({file_info.name})';

%% Loop through file to extract gains and lags for each lambda into a huge table
gainlag_cell = cell(length(filenames),1);
filetic = tic;
for filenum = 1:length(filenames)
    td = load(fullfile(dataroot,'library',filenames{filenum}));
    td = td.trial_data;

    [~,td_cst] = getTDidx(td,'task','CST','result','R');
    if ~isempty(td_cst)
        td_cst = getDifferential(td_cst,struct('signals','hand_pos','alias','hand_vel'));
        td_cst = getDifferential(td_cst,struct('signals','cursor_pos','alias','cursor_vel'));
    
        td_cst = trimTD(td_cst,{'idx_goCueTime',100},{'idx_rewardTime',-100});

        gainlag_td_cell = cell(length(td_cst),1);
        for trialnum = 1:length(td_cst)
            [corr_trace,corr_lag] = xcorr(td_cst(trialnum).hand_vel(:,1),td_cst(trialnum).cursor_vel(:,1));
            [~,lag_idx] = min(corr_trace);
            
            y_motion = sqrt(mean(td_cst(trialnum).hand_pos(:,2).^2));
            x_motion = sqrt(mean(td_cst(trialnum).hand_pos(:,1).^2));
            gain = sqrt(mean(td_cst(trialnum).hand_vel(:,1).^2))/...
                sqrt(mean(td_cst(trialnum).cursor_vel(:,1).^2));
            
            gainlag_td_cell{trialnum} = table(...
                {td_cst(trialnum).monkey},...
                {td_cst(trialnum).date_time},...
                td_cst(trialnum).trial_id,...
                td_cst(trialnum).lambda,...
                gain,...
                corr_lag(lag_idx)*td_cst(trialnum).bin_size,...
                x_motion,...
                y_motion,...
                'VariableNames',{'monkey','date_time','trial_id','lambda','gain','lag','x_motion','y_motion'});
            gainlag_td_cell{trialnum}.Properties.VariableDescriptions = {'meta','meta','meta','meta','linear','linear','linear','linear'};
        end
        gainlag_cell{filenum} = vertcat(gainlag_td_cell{:});
        fprintf('Finished file %d of %d at time %f\n',filenum,length(filenames),toc(filetic))
    else
        fprintf('No CST trials for file %d\n',filenum)
    end

end
gainlag_table = vertcat(gainlag_cell{:});

%% Select an intermediate lambda and plot progress
lambda_specific = 2.4;
[~,gainlag_specific] = getNTidx(gainlag_table,'monkey','Ford','lambda',lambda_specific);

% gainlag_specific(gainlag_specific.gain<1,:) = [];
% they should already be sorted, so just plot out lags
figure('defaultaxesfontsize',18)
subplot(2,1,1)
scatter(1:height(gainlag_specific),gainlag_specific.lag,[],'k','filled')
xlabel('Trial')
ylabel('Lag (s)')
set(gca,'box','off','tickdir','out')
subplot(2,1,2)
scatter(1:height(gainlag_specific),gainlag_specific.gain,[],'k','filled')
xlabel('Trial')
ylabel('Gain')
set(gca,'box','off','tickdir','out')
suptitle(sprintf('lambda = %0.1f',lambda_specific))

%% X and Y motion for intermediate lambda
lambda_specific = 3.2;
[~,gainlag_specific] = getNTidx(gainlag_table,'monkey','Ford','lambda',lambda_specific);

figure('defaultaxesfontsize',18)
subplot(2,1,1)
scatter(1:height(gainlag_specific),gainlag_specific.x_motion,[],'k','filled')
xlabel('Trial')
ylabel('X Motion RMS')
set(gca,'box','off','tickdir','out')
subplot(2,1,2)
scatter(1:height(gainlag_specific),gainlag_specific.y_motion,[],'k','filled')
xlabel('Trial')
ylabel('Y Motion RMS')
set(gca,'box','off','tickdir','out')
suptitle(sprintf('lambda = %0.1f',lambda_specific))

%% lambda vs lag
figure('defaultaxesfontsize',18)
scatter(gainlag_table.lambda,gainlag_table.lag,[],'k','filled')
xlabel('\lambda')
ylabel('Lag (s)')
set(gca,'box','off','tickdir','out')

% all trials
figure('defaultaxesfontsize',18)
subplot(2,1,1)
scatter(1:height(gainlag_table),gainlag_table.lag,[],gainlag_table.lambda,'filled')
colormap(viridis)
c = colorbar;
c.Label.String = '\lambda';
xlabel('Trial number')
ylabel('Lag (s)')
set(gca,'box','off','tickdir','out')
subplot(2,1,2)
scatter(1:height(gainlag_table),gainlag_table.gain,[],gainlag_table.lambda,'filled')
colormap(viridis)
c = colorbar;
c.Label.String = '\lambda';
xlabel('Trial number')
ylabel('Gain')
set(gca,'box','off','tickdir','out')

% 3d?
figure('defaultaxesfontsize',18)
scatter3(datenum(gainlag_table.date_time),gainlag_table.lambda,gainlag_table.lag,[],gainlag_table.trial_id,'filled')
colormap(viridis)
c = colorbar;
c.Label.String = 'Trial number';
xlabel('Date number')
ylabel('\lambda')
zlabel('Lag (s)')
set(gca,'box','off','tickdir','out')
