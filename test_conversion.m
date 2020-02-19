%% convert data
clearvars -except smile_data
td = convertSMILEtoTD(smile_data);

%% plot cursor positions for some Center Out trials
[~,td_co] = getTDidx(td,'task','CO','result','R');
td_co = trimTD(td_co,{'idx_goCueTime',0},{'idx_rewardTime',0});
ot_locs = unique(cat(1,td_co.ot_location),'rows');
dir_colors = linspecer(size(ot_locs,1));

figure
for trialnum = 1:length(td_co)
    dir_idx = ismember(ot_locs,td_co(trialnum).ot_location,'rows');
    plot(td_co(trialnum).hand_pos(:,1),td_co(trialnum).hand_pos(:,2),'color',dir_colors(dir_idx,:))
    hold on
end
axis equal
set(gca,'box','off','tickdir','out')

%% plot example CST trials
[~,td_cst] = getTDidx(td,'task','CST','result','R');
td_cst = trimTD(td_cst,{'idx_goCueTime',100},{'idx_rewardTime',-100});

figure
for trialnum = 1:30
    clf
    scatter(-td_cst(trialnum).ct_location(1),-td_cst(trialnum).ct_location(2),100,'k','filled')
    hold on
    plot(-td_cst(trialnum).hand_pos(:,1),-td_cst(trialnum).hand_pos(:,2),'-k','linewidth',2)
    set(gca,'box','off','tickdir','out')
    waitforbuttonpress
end

figure
for trialnum = 1:30
    clf
    plot(td_cst(trialnum).hand_pos(:,1),'-k','linewidth',2)
    hold on
    plot(td_cst(trialnum).cursor_pos(:,1),'-g','linewidth',2)
    set(gca,'box','off','tickdir','out')
    waitforbuttonpress
end

%% calculate cross-correlation between hand and cursor
[~,td_cst] = getTDidx(td,'task','CST','result','R');
td_cst = trimTD(td_cst,{'idx_goCueTime',100},{'idx_rewardTime',-100});

% figure
[gain,lag] = deal(zeros(length(td_cst),1));
for trialnum = 1:length(td_cst)
    [corr_trace,corr_lag] = xcorr(td_cst(trialnum).hand_pos(:,1),td_cst(trialnum).cursor_pos(:,1));
    [~,lag_idx] = min(corr_trace);
    lag(trialnum) = corr_lag(lag_idx)*td_cst(trialnum).bin_size;
    gain(trialnum) = corr_trace(lag_idx);

    % clf
    % plot(lags,corr_trace,'-k','linewidth',2)
    % set(gca,'box','off','tickdir','out')
    % waitforbuttonpress
end

figure
subplot(2,1,1)
scatter(1:length(td_cst),gain,[],'k','filled')
set(gca,'box','off','tickdir','out')
subplot(2,1,2)
scatter(1:length(td_cst),lag,[],'k','filled')
set(gca,'box','off','tickdir','out')
