%%    
    dataroot = '/data/raeed/project-data/smile/cst-gainlag';
    if ispc
        dataroot = 'C:\Users\Raeed\data\project-data\smile\cst-gainlag';
    end

    file_info = dir(fullfile(dataroot,'library','*COCST*'));
    filenames = horzcat({file_info.name})';
    
    gpfa_run_dir = 'C:\\Users\\Raeed\\data\\temp\\gpfa_runs';
    
%% load data
    file_query = struct(...
        'monkey','Ford',...
        'date','20180627');
    td_preproc = load_clean_cst_data(fullfile(dataroot,'library',sprintf('%s_%s_COCST_TD.mat',file_query.monkey,file_query.date)));
    
%% Smooth and extract td_cst
    td = td_preproc;
    [~,td_co] = getTDidx(td,'task','CO','result','R');
    [~,td_cst] = getTDidx(td,'task','CST','result','R');
    
    co_len = 500;
    cst_len = 5000;
    td_co = trimTD(td_co,{'idx_goCueTime',0},{'idx_goCueTime',co_len});
    td_cst = trimTD(td_cst,{'idx_cstStartTime',1000},{'idx_cstStartTime',1000+co_len});
    
%     num_cst_trials = floor(length(td_co)/(cst_len/co_len));
    num_trials = min(length(td_co),length(td_cst));
    [~,td_co] = getTDidx(td_co,'rand',num_trials);
    [~,td_cst] = getTDidx(td_cst,'rand',num_trials);
    
%     num_dims = 10;
    kernSD = 0.03;
    bin_w = 0.02;
    
%     td_binned = binTD(td,bin_w/td(1).bin_size);
    for num_dims = [2,5,8,10,15,20,25,30]
%         runGPFA(td_co,struct(...
%             'arrays','M1',...
%             'method','gpfa',...
%             'xDim',num_dims,...
%             'kernSD',kernSD,...
%             'bin_w',bin_w,...
%             'save_dir',gpfa_run_dir,...
%             'runid','-CO',...
%             'numFolds',5));
        
        runGPFA(td_cst,struct(...
            'arrays','M1',...
            'method','gpfa',...
            'xDim',num_dims,...
            'kernSD',kernSD,...
            'bin_w',bin_w,...
            'save_dir',gpfa_run_dir,...
            'runid','-CST-trialmatch',...
            'numFolds',5));
    end

%% show plot of dimensionality
    
    plotLLVsDim(gpfa_run_dir,'-M1-CO');
    set(gca,'box','off','tickdir','out')
    title('Center out LL vs. dims')
    plotLLVsDim(gpfa_run_dir,'-M1-CST');
    set(gca,'box','off','tickdir','out')
    title('CST LL vs. dims (full)')
    plotLLVsDim(gpfa_run_dir,'-M1-CST-trialmatch');
    set(gca,'box','off','tickdir','out')
    title('CST LL vs. dims (trial-matched)')
    plotLLVsDim(gpfa_run_dir,'-M1-CST-totmatch');
    set(gca,'box','off','tickdir','out')
    title('CST LL vs. dims (timepoint-matched)')
