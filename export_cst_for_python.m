%% Set up
    dataroot = '/data/raeed/project-data/smile/cst-gainlag';
    if ispc
        dataroot = 'C:\Users\Raeed\data\project-data\smile\cst-gainlag';
    end

    file_info = dir(fullfile(dataroot,'library','*COCST*'));
    filenames = horzcat({file_info.name})';
    
%% Loop through files
filetic = tic;
for filenum = 1%length(filenames)
    td = load_clean_cst_data(fullfile(dataroot,'library',filenames{filenum}));
    
    [~,td_cst]=getTDidx(td,'task','CST');
    
    % Make sure we have CST trials
    if isempty(td_cst)
        fprintf('Incomplete dataset for file %d\n',filenum)
        continue
    end
    
    % prep TD for phase plotting
%     td_cst = trimTD(td_cst,{'idx_cstStartTime',100},{'idx_cstEndTime',0});
    
    % add timevector to each trial
    for trialnum = 1:length(td_cst)
        td_cst(trialnum).trialtime = ((1:length(td_cst(trialnum).hand_pos))-1)'*td_cst(trialnum).bin_size;
    end
    
    save(fullfile(dataroot,'library','python',filenames{filenum}),'td_cst')

    fprintf('Finished file %d of %d at time %f\n',filenum,length(filenames),toc(filetic))
end