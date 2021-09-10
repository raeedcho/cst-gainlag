%%
    % dataroot = '/mnt/smile-share/Animals/';
    dataroot = 'C:\Users\Raeed\data\project-data\smile\cst-gainlag\';
    file_info = dir(fullfile(dataroot,'preTD','*CST*'));
    filenames = horzcat({file_info.name})';
    % savedir = '/data/raeed/project-data/smile/cst-gainlag/library';
    savedir = fullfile(dataroot,'library');

%%
%% go through files
filetic = tic;
for filenum = [1 3]
    smile_data = load(fullfile(dataroot,'preTD',filenames{filenum}));

    trial_data = convertSMILEtoTD(smile_data.Data,struct('array_alias',{{'Right M1','M1'}}));
    
    session_date = strsplit(trial_data(1).date_time);
    session_date = strrep(session_date{1},'/','');
    save(fullfile(savedir,sprintf('%s_%s_COCST_TD.mat',trial_data(1).monkey,session_date)),'trial_data','-v7.3')
    
    fprintf('Finished file %d of %d at time %f\n',filenum,length(filenames),toc(filetic));
end