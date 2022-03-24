% script to run Dekleva subspace splitter on CO and CST data

%% setup
    dataroot = '/data/raeed/project-data/smile/cst-gainlag';
    if ispc
        dataroot = 'C:\Users\Raeed\data\project-data\smile\cst-gainlag';
    end

    file_info = dir(fullfile(dataroot,'library','*COCST*'));
    filenames = horzcat({file_info.name})';
    
%% Load a file
    file_query = struct(...
        'monkey','Earl',...
        'date','20190716');
    td_preproc = load_clean_cst_data(fullfile(dataroot,'library',sprintf('%s_%s_COCST_TD.mat',file_query.monkey,file_query.date)));
    
    % Make sure we have CST trials
    assert(~isempty(td_preproc),sprintf('Incomplete dataset for file %s %s\n', file_query.monkey,file_query.date))
    assert(~isempty(td_preproc(1).M1_unit_guide),sprintf('Skipping file %s %s because no spike data...\n',file_query.monkey,file_query.date))

%% Preprocess td some more...
    % options
    analysis_epoch = 'move';
    smoothsigs = true;

    td = td_preproc;
    
    % smooth data
    if smoothsigs
        td = smoothSignals(td,struct('signals','M1_spikes','width',0.075,'calc_rate',true,'field_extra','_smooth'));
    end
    
    % split data
    [~,td_co] = getTDidx(td,'task','CO');
    [~,td_cst] = getTDidx(td,'task','CST');
    
    lambdas = sort(unique(cat(1,td_cst.lambda)));
    
    % trim TD to some arbitrary window after go cue
    if strcmpi(analysis_epoch,'move')
        td_co = trimTD(td_co,{'idx_goCueTime',150},{'idx_goCueTime',400});
        td_cst = trimTD(td_cst,{'idx_cstStartTime',150},{'idx_cstStartTime',5000});
    elseif strcmpi(analysis_epoch,'hold')
        td_co = trimTD(td_co,{'idx_goCueTime',-450},{'idx_goCueTime',0});
        td_cst = trimTD(td_cst,{'idx_goCueTime',-450},{'idx_goCueTime',0});
    end

    td_co = binTD(td_co,50);
    td_cst = binTD(td_cst,50);    

%% find joint space across both tasks
    num_dims = 10;

    % soft normalize and dim reduce separately
    td_co_norm = softNormalize(td_co,struct('signals','M1_spikes_smooth'));
    td_cst_norm = softNormalize(td_cst,struct('signals','M1_spikes_smooth'));
    [~,pca_co] = dimReduce(td_co_norm,struct(...
        'algorithm','pca',...
        'signals','M1_spikes_smooth',...
        'num_dims',num_dims));
    [~,pca_cst] = dimReduce(td_cst_norm,struct(...
        'algorithm','pca',...
        'signals','M1_spikes_smooth',...
        'num_dims',num_dims));

    % put spaces together and orthogonalize
    w_co = pca_co.w(:,1:num_dims);
    w_cst = pca_cst.w(:,1:num_dims);
    [w_new,~,~] = svd([w_co w_cst],0);

    % reproject data
    td_full = cat(2,td_co,td_cst);
    td_full = softNormalize(td_full,struct('signals','M1_spikes_smooth'));
    td_full = centerSignals(td_full,struct('signals','M1_spikes_smooth'));
    for trialnum = 1:length(td_full)
        td_full(trialnum).M1_pca_full = td_full(trialnum).M1_spikes_smooth * w_new;
    end

    % find private and shared dimensions of neural data across tasks
    % split data again
    [~,td_co] = getTDidx(td_full,'task','CO');
    [~,td_cst] = getTDidx(td_full,'task','CST');
    M1_state_co = getSig(td_co,'M1_pca_full');
    M1_state_cst = getSig(td_cst,'M1_pca_full');

    alpha_null_space = 1e-4;
    var_cutoff = 0.025;
    do_plot = true;

    [subspaces,projs,out] = SubspaceSplitterKaoOrthGreedy(M1_state_cst,M1_state_co,alpha_null_space,var_cutoff,do_plot);

    % relabel plots
    subplot(3,3,[1 2])
    legend({'CST','CO'},'Location','Best','box','off','FontSize',14); 
    subplot(3,3,3)
    xlabel('% var CST'); ylabel('% var CO'); 
    subplot(3,3,6)
    xlabel('% var CST'); ylabel('% var CO'); 
    subplot(3,3,[7 8])
    xticklabels({'CST unique','shared','CO unique'}); 
