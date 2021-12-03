% look at CST neural trajectories

%% Set up
    dataroot = '/data/raeed/project-data/smile/cst-gainlag';
    if ispc
        dataroot = 'C:\Users\Raeed\data\project-data\smile\cst-gainlag';
    end

    file_info = dir(fullfile(dataroot,'library','*COCST*'));
    filenames = horzcat({file_info.name})';
    
%% Select a file
    file_query = struct(...
        'monkey','Ford',...
        'date','20180618');
    td_preproc = load_clean_cst_data(fullfile(dataroot,'library',sprintf('%s_%s_COCST_TD.mat',file_query.monkey,file_query.date)));
    
    % Make sure we have CST trials
    assert(~isempty(td_preproc),sprintf('Incomplete dataset for file %s %s\n', file_query.monkey,file_query.date))
    assert(~isempty(td_preproc(1).M1_unit_guide),sprintf('Skipping file %s %s because no spike data...\n',file_query.monkey,file_query.date))
    
%%
    num_dims = 8;
    td = td_preproc;

    % trim TD to only center hold portion
    td = trimTD(td,{'idx_goCueTime',-450},{'idx_goCueTime',0});

    % smooth data
%     td = smoothSignals(td,struct('signals','M1_spikes','width',0.075,'calc_rate',true));
    td = softNormalize(td,struct('signals','M1_spikes','alpha',5));
    td = dimReduce(td,struct('algorithm','pca','signals','M1_spikes','num_dims',num_dims));
%     td = getDifferential(td,struct('signals','M1_pca','alias','M1_pca_diff'));

    td_binned = binTD(td,'average');
    
    % split data
    [~,td_co] = getTDidx(td_binned,'task','CO');
    [~,td_cst] = getTDidx(td_binned,'task','CST');
    
    % make plot of hold time activity
    figure
    subplot(1,2,1)
    scatter3(...
        get_vars(td_co,{'M1_pca',1}),...
        get_vars(td_co,{'M1_pca',2}),...
        get_vars(td_co,{'M1_pca',3}),...
        [],'r','filled');
    hold on
    scatter3(...
        get_vars(td_cst,{'M1_pca',1}),...
        get_vars(td_cst,{'M1_pca',2}),...
        get_vars(td_cst,{'M1_pca',3}),...
        [],cat(1,td_cst.lambda),'filled');
    colormap(viridis);
    set(gca,'box','off','tickdir','out')
    axis equal
    title('M1 PCA - CO in red, CST in viridis')
    
    subplot(1,2,2)
    scatter3(...
        get_vars(td_co,{'hand_pos',1}),...
        get_vars(td_co,{'hand_pos',2}),...
        get_vars(td_co,{'hand_pos',3}),...
        [],'r','filled');
    hold on
    scatter3(...
        get_vars(td_cst,{'hand_pos',1}),...
        get_vars(td_cst,{'hand_pos',2}),...
        get_vars(td_cst,{'hand_pos',3}),...
        [],cat(1,td_cst.lambda),'filled');
    set(gca,'box','off','tickdir','out')
    axis equal
    title('Hand position - CO in red, CST in viridis')
