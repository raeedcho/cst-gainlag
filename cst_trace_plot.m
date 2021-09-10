% look at CST/CO neural trajectories

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
start_time = -0.4;
end_time = 0.4;
filetic = tic;
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
    td = smoothSignals(td,struct('signals','M1_spikes','width',0.05,'calc_rate',true));
%     td = softNormalize(td,struct('signals','M1_spikes','alpha',5));
    td = dimReduce(td,struct('algorithm','pca','signals','M1_spikes','num_dims',num_dims));
    [td.M1_pca_full] = deal(td.M1_pca);
%     td = getDifferential(td,struct('signals','M1_pca','alias','M1_pca_diff'));

    % trim to time around go cue
    td = trimTD(td,struct(...
        'idx_start',{{'idx_goCueTime',floor(start_time/td(1).bin_size)}},...
        'idx_end',{{'idx_goCueTime',floor(end_time/td(1).bin_size)}},...
        'remove_short',true));
    
    % split data
    [~,td_co] = getTDidx(td,'task','CO');
    [~,td_cst] = getTDidx(td,'task','CST');
    
    % apply pca individually to CST and CO
    td_co = dimReduce(td_co,struct('algorithm','pca','signals','M1_spikes','num_dims',num_dims));
    td_cst = dimReduce(td_cst,struct('algorithm','pca','signals','M1_spikes','num_dims',num_dims));
    
    % plot out neural traces
    dirs = unique(cat(1,td_co.tgtDir));
    dir_colors = linspecer(length(dirs));
    figure
    for dirnum = 1:length(dirs)
        trial_idx = getTDidx(td_co,'task','CO','tgtDir',dirs(dirnum));
        plot_traces(td_co,struct(...
            'signals',{{'M1_pca_full',1:2}},...
            'trials_to_use',trial_idx,...
            'trials_to_plot',trial_idx(randperm(length(trial_idx),10)),...
            'plot_dim',2,...
            'color',dir_colors(dirnum,:)))
        hold on
    end
    trial_idx = getTDidx(td_cst,'task','CST','lambda',lambda_to_use,'trial_id',159);
    plot_traces(td_cst,struct(...
        'signals',{{'M1_pca_full',1:3}},...
        'trials_to_use',trial_idx,...
        'trials_to_plot',trial_idx(randperm(length(trial_idx),1)),...
        'plot_dim',2,...
        'color',[0 0 0]))
    
    set(gca,'box','off','tickdir','out')
    axis equal
    axis off
    
    % plot fake axes
%     xlims = get(gca,'xlim');
%     ylims = get(gca,'ylim');
%     zlims = get(gca,'zlim');
%     ax_len = diff(xlims)*0.2;
%     plot3([xlims(1) xlims(1)+ax_len],[ylims(1) ylims(1)],[zlims(1) zlims(1)],'k','linewidth',3)
%     plot3([xlims(1) xlims(1)],[ylims(1) ylims(1)+ax_len],[zlims(1) zlims(1)],'k','linewidth',3)
%     plot3([xlims(1) xlims(1)],[ylims(1) ylims(1)],[zlims(1) zlims(1)+ax_len],'k','linewidth',3)

    fprintf('Finished file %d of %d at time %f\n',filenum,length(filenames),toc(filetic))
end
