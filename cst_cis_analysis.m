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
    
    % smooth data
    td = smoothSignals(td,struct('signals','M1_spikes','width',0.075,'calc_rate',true));
%     td = softNormalize(td,struct('signals','M1_spikes','alpha',5));
    td = dimReduce(td,struct('algorithm','pca','signals','M1_spikes','num_dims',num_dims));
%     td = getDifferential(td,struct('signals','M1_pca','alias','M1_pca_diff'));
    
    % split data
    [~,td_co] = getTDidx(td,'task','CO');
    [~,td_cst] = getTDidx(td,'task','CST');
    
    % find condition independent signal for CO trials
    td_co = trimTD(td_co,{'idx_goCueTime',-450},{'idx_goCueTime',400});
    
    [td_co,dpca_info] = runDPCA(td_co,'tgtDir',struct(...
        'signals','M1_spikes',...
        'num_dims',[2 6],...
        'do_plot',false,...
        'marg_names',{{'time','target'}},...
        'out_sig_prefix','M1_dpca'));
%     co_fr = cat(3,td_co.M1_spikes);
%     co_fr = mean(co_fr,3);
%     co_fr = co_fr';
%     [W, V, which_marg] = dpca( co_fr, 2, 'combinedParams', {{1}} );
%     expl_var = dpca_explainedVariance(co_fr, W, V, ...
%         'combinedParams', {{1}});
    
    % make plot of CIS activity
    figure
    M1_time = cat(3,td_co.M1_dpca_time);
    plot(-0.45:0.001:0.4,mean(M1_time(:,2,:),3)')

    fprintf('Finished file %d of %d at time %f\n',filenum,length(filenames),toc(filetic))
end
