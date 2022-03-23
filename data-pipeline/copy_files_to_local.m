%%% move Ford files to project-data
%% Set up
dataroot = 'Y:\\Animals';
monkey = 'Ford';
fid = fopen('ford-cst-dates.txt');
dates = textscan(fid,'%s');
fclose(fid);
dates=dates{1};
% savedir = '/data/raeed/project-data/smile/cst-gainlag/library';
savedir = 'C:\\Users\\Raeed\\data\\project-data\\smile\\cst-gainlag\\preTD';

%% go through files
filetic = tic;
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

    copyfile(fullfile(file_info.folder,file_info.name),fullfile(savedir,sprintf('%s_%s_COCST_smile.mat',monkey,dates{filenum})))
    
    fprintf('Finished file %d of %d at time %f\n',filenum,length(dates),toc(filetic));
end