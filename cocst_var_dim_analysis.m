% look at CST neural trajectories

%% Set up
    dataroot = '/data/raeed/project-data/smile/cst-gainlag';
    if ispc
        dataroot = 'C:\Users\Raeed\data\project-data\smile\cst-gainlag';
    end

    file_info = dir(fullfile(dataroot,'library','*COCST*'));
    filenames = horzcat({file_info.name})';
    
%% Loop through files
filetic = tic;
num_dims = 8;
for filenum = 1%length(filenames)
    td = load_clean_cst_data(fullfile(dataroot,'library',filenames{filenum}));
    
    % Make sure we have CST trials
    if isempty(td)
        fprintf('Incomplete dataset for file %d\n',filenum)
        continue
    end
    
    if isempty(td(1).M1_unit_guide)
        fprintf('Skipping file %d because no spike data...\n',filenum)
        continue
    end
    
    % smooth data
    td = smoothSignals(td,struct('signals','M1_spikes','width',0.075,'calc_rate',true));
%     td = softNormalize(td,struct('signals','M1_spikes','alpha',5));
    td = dimReduce(td,struct('algorithm','pca','signals','M1_spikes','num_dims',num_dims));
%     td = getDifferential(td,struct('signals','M1_pca','alias','M1_pca_diff'));
    
    % split data
    [~,td_co] = getTDidx(td,'task','CO');
    [~,td_cst] = getTDidx(td,'task','CST');
    
    lambdas = sort(unique(cat(1,td_cst.lambda)));
    
    % trim TD to some arbitrary window after go cue
    td_co = trimTD(td_co,{'idx_goCueTime',150},{'idx_goCueTime',400});
    td_cst = trimTD(td_cst,{'idx_cstStartTime',150},{'idx_cstStartTime',400});
    
    % Get variances
    total_co_var = sum(var(getSig(td_co,'M1_spikes')));
    total_co_pr = calc_participation_ratio(getSig(td_co,'M1_spikes'));
    total_cst_var = sum(var(getSig(td_cst,'M1_spikes')));
    total_cst_pr = calc_participation_ratio(getSig(td_cst,'M1_spikes'));
    lambda_var = zeros(length(lambdas),1);
    lambda_pr = zeros(length(lambdas),1);
    for lambdanum = 1:length(lambdas)
        [~,td_lambda] = getTDidx(td_cst,'lambda',lambdas(lambdanum));
        lambda_var(lambdanum) = sum(var(getSig(td_lambda,'M1_spikes')));
        lambda_pr(lambdanum) = calc_participation_ratio(getSig(td_lambda,'M1_spikes'));
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
    set(gca,'box','off','tickdir','out','xlim',[min(lambdas) max(lambdas)],'ylim',[0 30])

    fprintf('Finished file %d of %d at time %f\n',filenum,length(filenames),toc(filetic))
end
