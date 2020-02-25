%% Set up
dataroot = '/mnt/smile-share/Animals/';
monkey = 'Ford';
fid = fopen('ford-cst-dates.txt');
dates = textscan(fid,'%s');
fclose(fid);
dates=dates{1};
savedir = '/data/raeed/project-data/smile/cst-gainlag/library';

%% go through files
for filenum = 1:length(dates)
    file_info = dir(fullfile(...
        dataroot,...
        monkey,...
        dates{filenum}(1:4),...
        dates{filenum}(5:6),...
        dates{filenum},...
        'handControl',...
        sprintf('%s*.mat',monkey)));

    if length(file_info)>1
        fprintf('Multiple matching files for date %s, using the first one, called %s...\n',dates{filenum},file_info(1).name)
        file_info = file_info(1);
    elseif isempty(file_info)
        fprintf('No file for date %s, skipping to next\n',dates{filenum})
        fid = fopen('ford-cst-untranslated.txt','a');
        fprintf(fid,'%s\n',dates{filenum});
        fclose(fid);
        continue
    end

    smile_data = load(fullfile(file_info.folder,file_info.name));

    trial_data = convertSMILEtoTD(smile_data.Data,struct('array_alias',{{'Right M1','M1'}}));
    save(fullfile(savedir,sprintf('%s_%s_TD.mat',monkey,dates{filenum})),'trial_data','-v7.3')
end
