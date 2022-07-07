%% plot raw analog signals out
% [~,td_co] = getTDidx(trial_data,'task','CO','result','R');
% [~,td_cst] = getTDidx(trial_data,'task','CST','result','R');
% td_co = trimTD(td_co,{'idx_goCueTime',-400},{'idx_endTime',0});
for trialnum=1:length(trial_data)
%     trialnum = 7;
    trial = trial_data(trialnum);
    h=figure;
    h_x_pos = plot(trial.raw_left_eye_pos(:,1));
    hold on
    h_y_pos = plot(trial.raw_left_eye_pos(:,2));
    h_pupil = plot(trial.raw_left_pupil);
    

    y_limits = get(gca,'ylim');
    plot([1;1]*trial.idx_ctHoldTime,y_limits,'--k','linewidth',2)
    plot([1;1]*trial.idx_goCueTime,y_limits,'--k','linewidth',2)
    plot([1;1]*trial.idx_rewardTime,y_limits,'--k','linewidth',2)
    plot([0 1000],[1 1]*y_limits(1)*0.5,'-k','linewidth',4)
    text(500,y_limits(1)*0.5-1,'1 second','FontSize',10)


    title(sprintf('Direction %.2f, Magnitude %.2f',trial.tgtDir,trial.tgtMag))
    ylabel('Channel voltage (V)')
    xlabel('Time in trial')
    legend([h_x_pos,h_y_pos,h_pupil],'Eye position X', 'Eye position Y', 'Pupil diameter')

    set(gca, ...
        'box','off', ...
        'tickdir','out', ...
        'xtick',[ ...
            trial.idx_ctHoldTime, ...
            trial.idx_goCueTime, ...
            trial.idx_rewardTime ...
            ],...
        'xticklabel',[ ...
            repmat({'Center Hold'},size(trial.idx_ctHoldTime)), ...
            repmat({'Go Cue'},size(trial.idx_goCueTime)), ...
            repmat({'Reward'},size(trial.idx_rewardTime)), ...
            ]...
        )
    
    waitfor(h)
end

%% See if gain scale makes sense and rescale eye data
% [idx_left,~] = getTDidx(trial_data,'task','CO','result','R','tgtDir',180);
% [idx_right,~] = getTDidx(trial_data,'task','CO','result','R','tgtDir',0);
% td_co = trial_data(union(idx_left,idx_right));
[~,td_co] = getTDidx(trial_data,'task','CO','result','R');
td_co = trimTD(td_co,{'idx_otHoldTime',-100},{'idx_otHoldTime',49});
eye_pos = cat(3,td_co.raw_left_eye_pos);
hand_pos = cat(3,td_co.hand_pos);
target_loc = repmat(cat(3,td_co.ot_location),[150,1,1]);
figure
scatter(reshape(hand_pos(:,2,:),[],1),reshape(eye_pos(:,2,:),[],1))
hold on
scatter(reshape(target_loc(:,2,:),[],1),reshape(eye_pos(:,2,:),[],1))

% coefs_x = polyfit(reshape(eye_pos(:,1,:),[],1),reshape(hand_pos(:,1,:),[],1),1);
% coefs_y = polyfit(reshape(eye_pos(:,2,:),[],1),reshape(hand_pos(:,2,:),[],1),1);
% From Earl 20190716 log file
coefs_x = [20 80];
coefs_y = [20 -415];

for trialnum = 1:length(trial_data)
    trial_data(trialnum).left_eye_pos = horzcat(...
        polyval(coefs_x,trial_data(trialnum).raw_left_eye_pos(:,1)),...
        polyval(coefs_y,trial_data(trialnum).raw_left_eye_pos(:,2))...
    );
end

[~,td_co] = getTDidx(trial_data,'task','CO','result','R');
td_co = trimTD(td_co,{'idx_otHoldTime',-100},{'idx_otHoldTime',49});
eye_pos = cat(3,td_co.left_eye_pos);    
hand_pos = cat(3,td_co.hand_pos);
target_loc = repmat(cat(3,td_co.ot_location),[150,1,1]);
figure
scatter(reshape(hand_pos(:,2,:),[],1),reshape(eye_pos(:,2,:),[],1))
hold on
scatter(reshape(target_loc(:,2,:),[],1),reshape(eye_pos(:,2,:),[],1))
plot(xlim,xlim,'--k')
xlabel('Screen position (hand or target)')
ylabel('Estimated eye position')

%% plot calibrated behavioral signals
% [~,td_co] = getTDidx(trial_data,'task','CO','result','R');
% [~,td_cst] = getTDidx(trial_data,'task','CST','result','R');
% td_co = trimTD(td_co,{'idx_goCueTime',-400},{'idx_endTime',0});
for trialnum=100:length(trial_data)
%     trialnum = 7;
    trial = trial_data(trialnum);
    h=figure;
    h_x_pos = plot(trial.left_eye_pos(:,1),'g');
    hold on
    h_y_pos = plot(trial.left_eye_pos(:,2),'g');
    h_hand_x_pos = plot(trial.hand_pos(:,1),'r');
    h_hand_y_pos = plot(trial.hand_pos(:,2),'r');
    h_cursor_x_pos = plot(trial.cursor_pos(:,1),'b');
    h_cursor_y_pos = plot(trial.cursor_pos(:,2),'b');
    

    y_limits = get(gca,'ylim');
    plot([1;1]*trial.idx_ctHoldTime,y_limits,'--k','linewidth',2)
    plot([1;1]*trial.idx_goCueTime,y_limits,'--k','linewidth',2)
    plot([1;1]*trial.idx_rewardTime,y_limits,'--k','linewidth',2)
    plot([0 1000],[1 1]*y_limits(1)*0.5,'-k','linewidth',4)
    text(500,y_limits(1)*0.5-1,'1 second','FontSize',10)


    title(sprintf('Direction %.2f, Magnitude %.2f',trial.tgtDir,trial.tgtMag))
    ylabel('Channel voltage (V)')
    xlabel('Time in trial')
    legend([h_x_pos,h_hand_x_pos,h_cursor_x_pos],'Eye position (X and Y)', 'Hand position (X and Y)', 'Cursor position (X and Y)')

    set(gca, ...
        'box','off', ...
        'tickdir','out', ...
        'xtick',[ ...
            trial.idx_ctHoldTime, ...
            trial.idx_goCueTime, ...
            trial.idx_rewardTime ...
            ],...
        'xticklabel',[ ...
            repmat({'Center Hold'},size(trial.idx_ctHoldTime)), ...
            repmat({'Go Cue'},size(trial.idx_goCueTime)), ...
            repmat({'Reward'},size(trial.idx_rewardTime)), ...
            ]...
        )
    
    waitfor(h)
end

%% Plot behavioral signals on screen coordinates
for trialnum=100:length(trial_data)
%     trialnum = 7;
    trial = trial_data(trialnum);
    h=figure;
    h_pos = plot(trial.left_eye_pos(:,1),trial.left_eye_pos(:,2),'g');
    hold on
    plot(trial.hand_pos(:,1),trial.hand_pos(:,2),'r')
    plot(trial.cursor_pos(:,1),trial.cursor_pos(:,2),'b')
    axis equal

    title(sprintf('Direction %.2f, Magnitude %.2f',trial.tgtDir,trial.tgtMag))
%     ylabel('Channel voltage (V)')
%     xlabel('Time in trial')
%     legend([h_x_pos,h_hand_x_pos,h_cursor_x_pos],'Eye position (X and Y)', 'Hand position (X and Y)', 'Cursor position (X and Y)')

    set(gca,'box','off','tickdir','out')
    
    waitfor(h)
end










