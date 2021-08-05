%%    
    dataroot = '/data/raeed/project-data/smile/cst-gainlag';
    if ispc
        dataroot = 'C:\Users\Raeed\data\project-data\smile\cst-gainlag';
    end

    file_info = dir(fullfile(dataroot,'library','*COCST*'));
    filenames = horzcat({file_info.name})';
    
%% load data
    filenum = 1;
    td_preproc = load_clean_cst_data(fullfile(dataroot,'library',filenames{filenum}));
    
%% Smooth and extract td_cst
    td = td_preproc;
    % smooth data
    td = smoothSignals(td,struct('signals','M1_spikes','width',0.075,'calc_rate',true));
    
%% Look for neural engagement signal in center hold period
    [~,td_cst] = getTDidx(td,'task','CST');
    td_hold = trimTD(td_cst,struct(...
        'idx_start',{{'idx_goCueTime',-425}},...
        'idx_end',{{'idx_goCueTime',0}},...
        'remove_short',true));
    td_hold = binTD(td_hold,'average');
    
    [td_hold,hold_pca_info] = dimReduce(td_hold,struct('algorithm','pca','signals','M1_spikes'));
    
    lambda_changes = find(diff(cat(1,td_hold.lambda))~=0)';
    failures = getTDidx(td_hold,'result','F');
    
    which_pc = 1;
    line_colors = linspecer(2);
    figure('defaultaxesfontsize',10)
    plot(repmat(failures,2,1)+0.5,[-25;25],'--','color',line_colors(2,:),'linewidth',1.5)
    hold on
    plot(repmat(lambda_changes+0.5,2,1),[-25;25],'--','color',line_colors(1,:),'linewidth',1.5)
    plot(1:length(td_hold),getSig(td_hold,{'M1_pca',which_pc}),'-k')
    scatter(1:length(td_hold),getSig(td_hold,{'M1_pca',which_pc}),[],cat(1,td_hold.lambda),'filled')
    xlabel('Trial number')
    ylabel(sprintf('Hold activity in PC %d',which_pc))
    set(gca,'box','off','tickdir','out')
    colormap(viridis)
    
    % Try to fit a linear model between some dimension of hold activity and the previous lambda
    td_hold(1).prev_lambda = td_hold(1).lambda;
    for trialnum = 2:length(td_hold)
        td_hold(trialnum).prev_lambda = td_hold(trialnum-1).lambda;
    end
    [td_hold,lambda_model] = getModel(td_hold,struct(...
        'model_type','linmodel',...
        'model_name','lambda_pred',...
        'in_signals',{{'M1_pca',1:16}},...
        'out_signals',{{'prev_lambda',1}}));
    
    point_color = repmat([0 0 0],length(td_hold),1);
    point_color(lambda_changes+2,:) = repmat([0 1 0],length(lambda_changes),1); % hold activity would change only after 
    figure('defaultaxesfontsize',10)
    scatter(cat(1,td_hold.prev_lambda),cat(1,td_hold.linmodel_lambda_pred),[],point_color,'filled')
    hold on
    plot([1.5 4.5],[1.5 4.5],'--k')
    xlabel('Actual previous \lambda')
    ylabel('Predicted previous \lambda (from hold activity)')
    set(gca,'box','off','tickdir','out')
    axis equal
    
%% look for behavioral subspace
    [~,td_co] = getTDidx(td,'task','CO');
    [~,td_cst] = getTDidx(td,'task','CST');
    
    td_cst = trimTD(td_cst,struct(...
        'idx_start',{{'idx_cstStartTime',100}},...
        'idx_end',{{'idx_cstEndTime',0}},...
        'remove_short',true));
    
    td_co = trimTD(td_co,struct(...
        'idx_start',{{'idx_goCueTime',100}},...
        'idx_end',{{'idx_goCueTime',500}},...
        'remove_short',true));
    
    [td_cst,cst_pca_info] = dimReduce(td_cst,struct('algorithm','pca','signals','M1_spikes'));
    [td_co,co_pca_info] = dimReduce(td_co,struct('algorithm','pca','signals','M1_spikes'));
    
    td_cst = binTD(td_cst,50);
    td_co = binTD(td_co,50);
    
    [td_cst,cst_hand_vel_decoder] = getModel(td_cst,struct(...
        'model_type','linmodel',...
        'model_name','cst_hand_vel_decoder',...
        'in_signals',{{'M1_pca',1:16}},...
        'out_signals',{{'hand_vel',1:2}}));
    
    [td_co,co_hand_vel_decoder] = getModel(td_co,struct(...
        'model_type','linmodel',...
        'model_name','co_hand_vel_decoder',...
        'in_signals',{{'M1_pca',1:16}},...
        'out_signals',{{'hand_vel',1:2}}));
    
    td_cst = getModel(td_cst,co_hand_vel_decoder);
    td_co = getModel(td_co,cst_hand_vel_decoder);
    
    task_to_plot={'co','cst'};
    figure
    for train_task = 1:length(task_to_plot)
        for test_task = 1:length(task_to_plot)
            subplot(2,2,(train_task-1)*length(task_to_plot)+test_task)
            eval(sprintf('test_td=td_%s',task_to_plot{test_task}))
            
            scatter(getSig(test_td,{'hand_vel',1}),getSig(test_td,{sprintf('linmodel_%s_hand_vel_decoder',task_to_plot{train_task}),1}),[],'k','filled')
            xlabel(sprintf('Actual hand velocity in %s',task_to_plot{test_task}))
            ylabel(sprintf('Predicted hand velocity by %s model',task_to_plot{train_task}))
            axis equal
            set(gca,'box','off','tickdir','out','xlim',[-250 350],'ylim',[-250 350])
        end
    end
    