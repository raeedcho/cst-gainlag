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
for filenum = 27%length(filenames)
    td_cst = load_clean_cst_data(fullfile(dataroot,'library',filenames{filenum}));
    
    % get rid of unsorted neurons
    bad_units = td_cst(1).M1_unit_guide(:,2)<=1;
    % for file 31 only! remove neurons 8,2 and 64,2
    if contains(td_cst(1).date_time,'2018/06/26')
        corr_units = [8 2;64 2];
        bad_units = bad_units | ismember(td_cst(1).M1_unit_guide,corr_units,'rows');
    end
    for trialnum = 1:length(td_cst)
        td_cst(trialnum).M1_spikes = td_cst(trialnum).M1_spikes(:,~bad_units);
        td_cst(trialnum).M1_unit_guide = td_cst(trialnum).M1_unit_guide(~bad_units,:);
    end
    
    if isempty(td_cst(1).M1_unit_guide)
        fprintf('Skipping file %d because no spike data...\n',filenum)
        continue
    end
    
    % Make sure we have CST trials
    if isempty(td_cst)
        fprintf('Incomplete dataset for file %d\n',filenum)
        continue
    end
    
    % smooth data
    td_cst = smoothSignals(td_cst,struct('signals','M1_spikes','width',0.075,'calc_rate',true));
    
    % trim TD to only CST portion
    td_cst = trimTD(td_cst,'idx_cstStartTime','idx_cstEndTime');
    
    % check tangling over all trials in a particular lambda
    lambda_to_use = 3.3;
    [~,td_lambda] = getTDidx(td_cst,'lambda',lambda_to_use);
    tangle_struct = 

    fprintf('Finished file %d of %d at time %f\n',filenum,length(filenames),toc(filetic))
end