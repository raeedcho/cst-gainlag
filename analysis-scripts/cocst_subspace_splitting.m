% script to run Dekleva subspace splitter on CO and CST data
function cocst_subspace_splitting(file_query,params)
%% setup
    dataroot = '/data/raeed/project-data/smile/cst-gainlag';
    if ispc
        dataroot = 'C:\Users\Raeed\data\project-data\smile\cst-gainlag';
    end
    
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
    num_dims = 10;
    equalize_samples = false;
    num_samples = 1500;
    avg_co_data = false;
%     assignParams(who,params);

    td = td_preproc;
    
    % smooth data
    if smoothsigs
        td = smoothSignals(td,struct('signals','M1_spikes','width',0.075,'calc_rate',true,'field_extra','_smooth'));
    end
    
    % split data
    [~,td_co] = getTDidx(td,'task','CO');
    [~,td_cst] = getTDidx(td,'task','CST');
    
    % trim TD to some arbitrary window after go cue
    if strcmpi(analysis_epoch,'move')
        td_co = trimTD(td_co,struct('idx_start',{{'idx_goCueTime',150}},'idx_end',{{'idx_goCueTime',500}},'remove_short',true));
        td_cst = trimTD(td_cst,struct('idx_start',{{'idx_cstStartTime',150}},'idx_end',{{'idx_cstStartTime',5000}},'remove_short',true));
    elseif strcmpi(analysis_epoch,'hold')
        td_co = trimTD(td_co,struct('idx_start',{{'idx_goCueTime',-450}},'idx_end',{{'idx_goCueTime',0}},'remove_short',true));
        td_cst = trimTD(td_cst,struct('idx_start',{{'idx_goCueTime',-450}},'idx_end',{{'idx_goCueTime',0}},'remove_short',true));
    end

    td_co = binTD(td_co,50);
    td_cst = binTD(td_cst,50);

    % get average CO activity
    td_co_avg = trialAverage(td_co,'tgtDir');
    deleted_fields = setdiff(fieldnames(td_cst),fieldnames(td_co_avg));
    for fn = 1:length(deleted_fields)
        for trialnum = 1:length(td_co_avg)
            td_co_avg(trialnum).(deleted_fields{fn}) = NaN;
        end
    end
    if avg_co_data
        td_co = td_co_avg;
    end

    [td_cell,~] = joint_dim_reduce({td_cst,td_co},struct('signals','M1_spikes_smooth','num_dims',num_dims,'combine_before_projection',false));
    [td_cst,td_co] = deal(td_cell{:});

    % test...
    assert(length(unique({td_co.task}))==1)
    assert(length(unique({td_cst.task}))==1)
    
    M1_state_cst = getSig(td_cst,'M1_pca_joint');
    M1_state_co = getSig(td_co,'M1_pca_joint');

    % take random samples of each set of states to make things even
    if equalize_samples
        M1_state_cst = datasample(M1_state_cst,num_samples,'Replace',true);
        M1_state_co = datasample(M1_state_co,num_samples,'Replace',true);
    end

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

    % plot out dynamics in each subspace for CST trials? % variance in each subspace as a function of time?
    
end


