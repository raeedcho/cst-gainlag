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
    
%     % get rid of unsorted neurons
%     bad_units = td(1).M1_unit_guide(:,2)<=1;
%     % for file 31 only! remove neurons 8,2 and 64,2
%     if contains(td(1).date_time,'2018/06/26')
%         corr_units = [8 2;64 2];
%         bad_units = bad_units | ismember(td(1).M1_unit_guide,corr_units,'rows');
%     end
%     for trialnum = 1:length(td)
%         td(trialnum).M1_spikes = td(trialnum).M1_spikes(:,~bad_units);
%         td(trialnum).M1_unit_guide = td(trialnum).M1_unit_guide(~bad_units,:);
%     end
%     
%     if isempty(td(1).M1_unit_guide)
%         fprintf('Skipping file %d because no spike data...\n',filenum)
%         continue
%     end
%     
    % Make sure we have CST trials
    if isempty(td_cst)
        fprintf('Incomplete dataset for file %d\n',filenum)
        continue
    end
    
    % prep TD for phase plotting
    td_cst = trimTD(td_cst,'idx_cstStartTime','idx_cstEndTime');
    td_cst = calcTolInstab(td_cst);
    td_cst = findRestorationBlocks(td_cst);

    % make interactive CST phase plot
    h = plot_interactive_cst_phase(td_cst,struct(...
        'cursor_sig',{{'cursor_pos',1}},...
        'hand_sig',{{'gainlag_model',1}}));
    close(h)

    fprintf('Finished file %d of %d at time %f\n',filenum,length(filenames),toc(filetic))
end