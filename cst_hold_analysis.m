% look at CST neural trajectories

%% Set up
    dataroot = '/data/raeed/project-data/smile/cst-gainlag';
    if ispc
        dataroot = 'C:\Users\Raeed\data\project-data\smile\cst-gainlag';
    end

    file_info = dir(fullfile(dataroot,'library','*COCST*'));
    filenames = horzcat({file_info.name})';
    
%% Loop through files
lambda_to_use = 3.3;
num_dims = 8;
filetic = tic;
for filenum = 27%length(filenames)
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
    
    % get rid of unsorted neurons
    bad_units = td(1).M1_unit_guide(:,2)<=1;
    % for file 31 only! remove neurons 8,2 and 64,2
    if contains(td(1).date_time,'2018/06/26')
        corr_units = [8 2;64 2];
        bad_units = bad_units | ismember(td(1).M1_unit_guide,corr_units,'rows');
    end
    for trialnum = 1:length(td)
        td(trialnum).M1_spikes = td(trialnum).M1_spikes(:,~bad_units);
        td(trialnum).M1_unit_guide = td(trialnum).M1_unit_guide(~bad_units,:);
    end
    
    td = removeBadNeurons(td,struct(...
        'do_fr_check',true,...
        'min_fr',0.1,...
        'fr_window',{{'idx_cstStartTime',0;'idx_cstEndTime',0}},...
        'calc_fr',true,...
        'use_trials',getTDidx(td,'lambda',lambda_to_use)));
    
    % smooth data
    td = smoothSignals(td,struct('signals','M1_spikes','width',0.075,'calc_rate',true));
%     td = softNormalize(td,struct('signals','M1_spikes','alpha',5));
    td = dimReduce(td,struct('algorithm','pca','signals','M1_spikes','num_dims',num_dims));
%     td = getDifferential(td,struct('signals','M1_pca','alias','M1_pca_diff'));

    % trim TD to only center hold portion
    td = trimTD(td,{'idx_goCueTime',-450},{'idx_goCueTime',0});
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

    fprintf('Finished file %d of %d at time %f\n',filenum,length(filenames),toc(filetic))
end
