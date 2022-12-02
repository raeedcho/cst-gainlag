%%
    % dataroot = '/mnt/smile-share/Animals/';
    dataroot = 'C:\Users\Raeed\data\project-data\smile\cst-gainlag\';
    file_info = dir(fullfile(dataroot,'preTD','Prez2022072*.mat'));
    filenames = horzcat({file_info.name})';
    % savedir = '/data/raeed/project-data/smile/cst-gainlag/library';
    savedir = fullfile(dataroot,'library');

%%
%% go through files
filetic = tic;
for filenum = 1:length(filenames)
    smile_data = load(fullfile(dataroot,'preTD',filenames{filenum}));

    trial_data = convertSMILEtoTD_prez(smile_data.Data);

    % only RTT and CST trials...
    trial_data = trial_data(strcmpi({trial_data.task},'RTT') | strcmpi({trial_data.task},'CST'));
    
    session_date = strsplit(trial_data(1).date_time);
    session_date = strrep(session_date{1},'/','');
%     save(fullfile(savedir,sprintf('%s_%s_RTTCST_TD.mat',trial_data(1).monkey,session_date)),'trial_data','-v7.3')
    
    fprintf('Finished file %d of %d at time %f\n',filenum,length(filenames),toc(filetic));
end