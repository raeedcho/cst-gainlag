%% Set up
    dataroot = '/data/raeed/project-data/smile/cst-gainlag';

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
        td_cst = trimTD(td_cst,{'idx_goCueTime',100},{'idx_rewardTime',-100});

        gainlag_td_cell = cell(length(td_cst),1);
        for trialnum = 1:length(td_cst)
            [corr_trace,corr_lag] = xcorr(td_cst(trialnum).hand_pos(:,1),td_cst(trialnum).cursor_pos(:,1));
            [~,lag_idx] = min(corr_trace);

            gainlag_td_cell{trialnum} = table(...
                {td_cst(trialnum).monkey},...
                {td_cst(trialnum).date_time},...
                td_cst(trialnum).trial_id,...
                td_cst(trialnum).lambda,...
                corr_trace(lag_idx),...
                corr_lag(lag_idx)*td_cst(trialnum).bin_size,...
                'VariableNames',{'monkey','date_time','trial_id','lambda','gain','lag'});
            gainlag_td_cell.Properties.VariableDescriptions = {'meta','meta','meta','meta','linear','linear'};
        end
        gainlag_cell{filenum} = vertcat(gainlag_td_cell{:});
        fprintf('Finished file %d at time %f\n',filenum,toc(filetic))
    else
        fprintf('No CST trials for file %d\n',filenum)
    end

end
gainlag_table = vertcat(gainlag_cell{:});

%% Select an intermediate lambda and plot progress
[~,gainlag_specific] = getNTidx(gainlag_table,'monkey','Ford','lambda',2.4);

% they should already be sorted, so just plot out lags
figure('defaultaxesfontsize',18)
scatter(1:height(gainlag_specific),gainlag_specific.lag,[],'k','filled')
xlabel('Trial')
ylabel('Lag (s)')
set(gca,'box','off','tickdir','out')

% lambda vs lag
figure('defaultaxesfontsize',18)
scatter(gainlag_table.lambda,gainlag_table.lag,[],'k','filled')
xlabel('\lambda')
ylabel('Lag (s)')
set(gca,'box','off','tickdir','out')

% all trials
figure('defaultaxesfontsize',18)
scatter(1:height(gainlag_table),gainlag_table.lag,[],gainlag_table.lambda,'filled')
colormap(viridis)
colorbar on
xlabel('Trial number')
ylabel('Lag (s)')
set(gca,'box','off','tickdir','out')
