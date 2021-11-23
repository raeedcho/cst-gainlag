% look at CST neural trajectories

%% Set up
    dataroot = '/data/raeed/project-data/smile/cst-gainlag';
    if ispc
        dataroot = 'C:\Users\Raeed\data\project-data\smile\cst-gainlag';
    end

    file_info = dir(fullfile(dataroot,'library','*COCST*'));
    filenames = horzcat({file_info.name})';
    
%% Loop through files
    file_query = struct(...
        'monkey','Earl',...
        'date','20190716');
    td_preproc = load_clean_cst_data(fullfile(dataroot,'library',sprintf('%s_%s_COCST_TD.mat',file_query.monkey,file_query.date)));
    
    % Make sure we have CST trials
    assert(~isempty(td_preproc),sprintf('Incomplete dataset for file %s %s\n', file_query.monkey,file_query.date))
    assert(~isempty(td_preproc(1).M1_unit_guide),sprintf('Skipping file %s %s because no spike data...\n',file_query.monkey,file_query.date))
    
%% 
    td = td_preproc;
    
    % smooth data
%     td = smoothSignals(td,struct('signals','M1_spikes','width',0.075,'calc_rate',true));
    
    % split data
    [~,td_co] = getTDidx(td,'task','CO');
    [~,td_cst] = getTDidx(td,'task','CST');
    
    lambdas = sort(unique(cat(1,td_cst.lambda)));
    
    % trim TD to some arbitrary window after go cue
    td_co = trimTD(td_co,{'idx_goCueTime',150},{'idx_goCueTime',400});
    td_cst = trimTD(td_cst,{'idx_cstStartTime',150},{'idx_cstStartTime',5000});
%     td_co = trimTD(td_co,{'idx_goCueTime',-450},{'idx_goCueTime',0});
%     td_cst = trimTD(td_cst,{'idx_goCueTime',-450},{'idx_goCueTime',0});

    td_co = binTD(td_co,50);
    td_cst = binTD(td_cst,50);
    
    % Get variances
    total_co_var = sum(var(getSig(td_co,'M1_spikes')));
    total_co_pr = calc_participation_ratio(getSig(td_co,'M1_spikes'));
    total_co_pa_dim = parallel_analysis_dimensionality(getSig(td_co,'M1_spikes'),struct('do_plot','true'));
    total_cst_var = sum(var(getSig(td_cst,'M1_spikes')));
    total_cst_pr = calc_participation_ratio(getSig(td_cst,'M1_spikes'));
    total_cst_pa_dim = parallel_analysis_dimensionality(getSig(td_cst,'M1_spikes'),struct('do_plot','true'));
    lambda_var = zeros(length(lambdas),1);
    lambda_pr = zeros(length(lambdas),1);
    lambda_pa_dim = zeros(length(lambdas),1);
    for lambdanum = 1:length(lambdas)
        [~,td_lambda] = getTDidx(td_cst,'lambda',lambdas(lambdanum));
        lambda_var(lambdanum) = sum(var(getSig(td_lambda,'M1_spikes')));
        lambda_pr(lambdanum) = calc_participation_ratio(getSig(td_lambda,'M1_spikes'));
        lambda_pa_dim(lambdanum) = parallel_analysis_dimensionality(getSig(td_lambda,'M1_spikes'));
    end
    
    % plot vars
    figure('defaultaxesfontsize',18)
    plot(lambdas,lambda_var,'-k','linewidth',2)
    hold on
    plot(lambdas,repmat(total_co_var,size(lambdas)),'--r','linewidth',2)
    plot(lambdas,repmat(total_cst_var,size(lambdas)),'--k','linewidth',2)
    xlabel('\lambda')
    ylabel('Neural variance')
    set(gca,'box','off','tickdir','out','xlim',[min(lambdas) max(lambdas)],'ylim',[0 8e3])
    
    % plot PRs
    figure('defaultaxesfontsize',18)
    plot(lambdas,lambda_pr,'-k','linewidth',2)
    hold on
    plot(lambdas,repmat(total_co_pr,size(lambdas)),'--r','linewidth',2)
    plot(lambdas,repmat(total_cst_pr,size(lambdas)),'--k','linewidth',2)
    xlabel('\lambda')
    ylabel('Dimensionality (via PR)')
    set(gca,'box','off','tickdir','out','xlim',[min(lambdas) max(lambdas)],'ylim',[0 35])
    
    % plot PA dims
    figure('defaultaxesfontsize',18)
    plot(lambdas,lambda_pa_dim,'-k','linewidth',2)
    hold on
    plot(lambdas,repmat(total_co_pa_dim,size(lambdas)),'--r','linewidth',2)
    plot(lambdas,repmat(total_cst_pa_dim,size(lambdas)),'--k','linewidth',2)
    xlabel('\lambda')
    ylabel('Dimensionality (via parallel ananlysis)')
    set(gca,'box','off','tickdir','out','xlim',[min(lambdas) max(lambdas)],'ylim',[0 35])
